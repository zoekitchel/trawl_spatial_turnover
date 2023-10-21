
###########################
#This scripts explores why our patterns for Scotland Sea do not match Magurran et al. 2015. Turns out, DATRAS updated data soon after publcation of Magurran et al, and that has led to the changes.
#######################
##VERSIONS##
#######################
#R 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
#macOS Big Sur 11.7

#######################
##PACKAGES##
#######################

library(data.table)
library(vegan)
library(sf)
library(concaveman) #polygon around points
library(betapart) #allows us to partition beta diversity
library(geosphere)
library(ggpubr) #stat_regline_equation
library(nlme)
library(mapplots)
library(betapart)
library(icesDatras)

#######################
##DATA##
#######################

#Pull in cleaned FishGlob data

FishGlob_clean <- readRDS(
  here::here("data", "cleaned", "FishGlob.10year.spp_manualclean.rds")) #EXCLUDES


#######################
##FURTHER CLEAN DATA AND REDUCE TO SCOTTISH SURVEYS
#######################


#delete all rows where wgt_cpue is NA,  unless it is MEDITS region
FishGlob_clean.noNA <- FishGlob_clean[(survey != "MEDITS" & !is.na(wgt_cpue)) | survey == "MEDITS"]

#check
FishGlob_clean.noNA <- copy(FishGlob_clean.noNA)[!is.infinite(wgt_cpue), ]


#Some spp are duplicated within a haul. In all cases,  this is due to verbatim names matching to a single accepted name. Therefore,  we will sum these wgts and abundances

FishGlob_clean.noNA[, wgt := sum(wgt,  na.rm = T), .(accepted_name,  haul_id)][, 
                   wgt_h := sum(wgt_h,  na.rm = T), .(accepted_name,  haul_id)][, 
                   wgt_cpue := sum(wgt_cpue,  na.rm = T), .(accepted_name,  haul_id)][, 
                   num := sum(num,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                                                                                        num_h := sum(num_h,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                                                                                                                                                     num_cpue := sum(num_cpue,  na.rm = T), .(accepted_name,  haul_id)]

FishGlob_clean.noNA <- unique(FishGlob_clean.noNA[,
     .(survey, haul_id, country, sub_area, continent, stat_rec, station, stratum, year, month,
       day, quarter, season, latitude, longitude, haul_dur, area_swept, gear, depth, sbt, 
       sst, num, num_h, num_cpue, wgt, wgt_h, wgt_cpue, verbatim_aphia_id, accepted_name, 
       aphia_id, SpecCode, kingdom, phylum, class, order, family, genus, rank, survey_unit)])

#limit to Scotland only                                                                                                                                                                                                                                                                       num_h := sum(num_h,  na.rm = T), .(accepted_name,  haul_id)][, 
scotland_shelf_sea <- FishGlob_clean.noNA[survey=="SWC-IBTS",]

rm(FishGlob_clean)
rm(FishGlob_clean.noNA)


###################
#CLEANING TO MATCH MAGURRAN 2015
##################
#We focus on the 35 rectangles (Supplementary Fig. 1) for which there are good community time series.
ICES_rec_include <- c(444,441,440,442,451,453,450,452,
454,462,463,473,474,475,464,465,
395,394,390,391,392,400,401,402,403,
431,433,422,420,421,413,423,412,410,411)

#but, for some reason this paper uses alternate naming system for ICES rectangles
#matched using this visual
#UK map
#create data frame of east coast CA and USA from rnaturalearth
uk_basemap <- rnaturalearth::ne_countries(country = "united kingdom", scale = "medium")

uk_basemap <- raster::crop(uk_basemap, extent(-10,-3, 54.5,70))

#first assign the rectangles to nine 0.5˚ latitudinal bands
xlim <- c(-10,-3)
ylim <- c(54.5,70)
basemap(xlim, ylim)
sp::plot(uk_basemap)
draw.rect()

#see screenshot "magurran_2015_ices_rec_match.png"

#
ICES_rec_include <- c("47E3","47E4","47E5",
                      "46E2","46E3","46E4","46E5",
                      "45E0","45E1","45E2","45E3","45E4",
                      "44E0","44E1","44E3","44E4",
                      "43E1","43E3",
                      "42E0","42E1","42E2","42E3",
                      "41E0","41E1","41E2","41E3",
                      "40E0","40E1","40E2","40E3",
                      "39E0","39E1","39E2","39E4","39E5")
                      

#As 1992 and 1995 have substantially fewer samples than other years, we exclude them from our analyses.
#check these are already excluded
unique(scotland_shelf_sea[,.(year)])
#yes

#Sample rarefaction (as in ref. 6) ensures equal sampling effort across bands and is used in the calculation of 
#temporal a diversity (species richness) and temporal b diversity (Jaccard similarity 
#(presence/absence) and Bray–Curtis similarity (quantitative)).
#####-> we will not do this, because it's not used in spatial beta diversity calculations

#split into latitudinal bands
scotland_shelf_sea[,stat_rec_lat := substr(stat_rec,1,2)]

#sum across latitudinal bands
scotland_shelf_sea[,num_cpue_latband :=
                     sum(num_cpue, na.rm = T),.(stat_rec_lat,year,survey_unit)][,
                    wgt_cpue_latband :=
                      sum(wgt_cpue, na.rm = T),.(stat_rec_lat,year,survey_unit)]


#summarize
scotland_shelf_sea.r <- unique(scotland_shelf_sea[,.(survey,stat_rec_lat,year,quarter,
                                                     wgt_cpue_latband, num_cpue_latband,
                                                     verbatim_aphia_id,accepted_name,
                                                     aphia_id, SpecCode, kingdom,
                                                     phylum, class, order, family,
                                                     genus,  rank, survey_unit)])


###############
#Calculate Beta Dissimilarity between all latitudes for Q1 and Q2 surveys
###############
#loop through all regions and years

distances_dissimilarities_scotland_shelf_wgt <- data.table("stat_rec_1" = integer(), 
                                                 "stat_rec_2" = integer(), 
                                                 "bray_curtis_dissimilarity_balanced" = numeric(), 
                                                 year = integer(), 
                                                 abund_biomass  =  character(), 
                                                 "jaccard_dissimilarity_turnover" = numeric(), 
                                                 "jaccard_dissimilarity_nestedness" =  numeric(), 
                                                 "bray_curtis_dissimilarity_gradient" = as.numeric(), 
                                                 "jaccard_dissimilarity_total" = numeric(), 
                                                 "bray_curtis_dissimilarity_total" = numeric(), 
                                                 survey = character(), 
                                                 survey_unit = character())

distances_dissimilarities_scotland_shelf_abun <- data.table("stat_rec_1" = integer(), 
                                                           "stat_rec_2" = integer(), 
                                                           "bray_curtis_dissimilarity_balanced" = numeric(), 
                                                           year = integer(), 
                                                           abund_biomass  =  character(), 
                                                           "jaccard_dissimilarity_turnover" = numeric(), 
                                                           "jaccard_dissimilarity_nestedness" =  numeric(), 
                                                           "bray_curtis_dissimilarity_gradient" = as.numeric(), 
                                                           "jaccard_dissimilarity_total" = numeric(), 
                                                           "bray_curtis_dissimilarity_total" = numeric(), 
                                                           survey = character(), 
                                                           survey_unit = character())

scotland_survey_units <- unique(scotland_shelf_sea[,factor(survey_unit)])

for (i in 1:length(scotland_survey_units)) {
  
  scotland_shelf_sea.r_subset <- scotland_shelf_sea.r[survey_unit == scotland_survey_units[i], ]
  
  years <- sort(unique(scotland_shelf_sea.r_subset[, as.numeric(year)]))
  
  #haul id keys
  lat_recs <- scotland_shelf_sea.r_subset[,.(stat_rec_lat)]
  
  for (j in 1:length(years)) {
    reduced_year <- scotland_shelf_sea.r_subset[year == years[j], ]
      
    
    #FIRST DO BY WGT
      #if some rows have wgt_cpue missing,  get rid of these rows
      reduced_year_wgt <- reduced_year[complete.cases(reduced_year[, wgt_cpue_latband]), ]
      
      reduced_year_wide <- dcast(reduced_year_wgt,  stat_rec_lat + year ~ accepted_name,  value.var = "wgt_cpue_latband",  fun.aggregate = sum) #longitude to wide data for community matrix,  column names are latitudinal band then species
      
      
      ncols <- ncol(reduced_year_wide)
      communitymatrix <- reduced_year_wide[, 3:ncols] #community matrix
      communitymatrix.occurence <- communitymatrix
      communitymatrix.occurence[communitymatrix.occurence > 0] <- 1
      
      #list of stat_rec_lats
      key_IDs_subset <- reduced_year_wide[,stat_rec_lat]
      
      dissimilarities_abundance <- beta.pair.abund(communitymatrix,  index.family = "bray") #dissimilarity 
      dissimilarities_occurrence <- beta.pair(communitymatrix.occurence,  index.family = "jaccard") #dissimilarity
      
      #make into matrix
      dissimilarities_abundance_balanced.m <- as.matrix(dissimilarities_abundance$beta.bray.bal,  labels=TRUE) #bal = balanced
      dissimilarities_abundance_gradient.m <- as.matrix(dissimilarities_abundance$beta.bray.gra,  labels=TRUE) #gra = gradient
      dissimilarities_abundance_total.m <- as.matrix(dissimilarities_abundance$beta.bray,  labels=TRUE) #total
      
      dissimilarities_occurrence_turnover.m <- as.matrix(dissimilarities_occurrence$beta.jtu,  labels=TRUE) #jtu = turnover
      dissimilarities_occurrence_nestedness.m <- as.matrix(dissimilarities_occurrence$beta.jne,  labels=TRUE) #jne = nestedness
      dissimilarities_occurrence_total.m <- as.matrix(dissimilarities_occurrence$beta.jac,  labels=TRUE) #total
      
      colnames(dissimilarities_abundance_balanced.m) <- rownames(dissimilarities_abundance_balanced.m) <- key_IDs_subset
      colnames(dissimilarities_abundance_gradient.m) <- rownames(dissimilarities_abundance_gradient.m) <- key_IDs_subset
      colnames(dissimilarities_abundance_total.m) <- rownames(dissimilarities_abundance_total.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_turnover.m) <- rownames(dissimilarities_occurrence_turnover.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_nestedness.m) <- rownames(dissimilarities_occurrence_nestedness.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_total.m) <- rownames(dissimilarities_occurrence_total.m) <- key_IDs_subset
      
      #reshape dissimilarities
      dissimilarities_abundance_balanced.l <- reshape2::melt(dissimilarities_abundance_balanced.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_balanced")
      dissimilarities_abundance_gradient.l <- reshape2::melt(dissimilarities_abundance_gradient.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_gradient")
      dissimilarities_abundance_total.l <- reshape2::melt(dissimilarities_abundance_total.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_total")
      
      dissimilarities_occurrence_turnover.l <- reshape2::melt(dissimilarities_occurrence_turnover.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_turnover")
      dissimilarities_occurrence_nestedness.l <- reshape2::melt(dissimilarities_occurrence_nestedness.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_nestedness")
      dissimilarities_occurrence_total.l <- reshape2::melt(dissimilarities_occurrence_total.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_total")
      
      #and then to data table format
      dissimilarities_abundance_balanced.l <- data.table(dissimilarities_abundance_balanced.l)
      dissimilarities_abundance_gradient.l <- data.table(dissimilarities_abundance_gradient.l)
      dissimilarities_abundance_total.l <- data.table(dissimilarities_abundance_total.l)
      dissimilarities_occurrence_turnover.l <- data.table(dissimilarities_occurrence_turnover.l)
      dissimilarities_occurrence_nestedness.l <- data.table(dissimilarities_occurrence_nestedness.l)
      dissimilarities_occurrence_total.l <- data.table(dissimilarities_occurrence_total.l)
      
      #add year for these values
      dissimilarities_abundance_balanced.l[,  "year" := years[j]]
      dissimilarities_occurrence_turnover.l[,  "year" := years[j]]
      
      
      dissimilarities_abundance_gradient.l[,  "year" := years[j]]
      dissimilarities_occurrence_nestedness.l[,  "year" := years[j]]
      
      dissimilarities_abundance_total.l[,  "year" := years[j]]
      dissimilarities_occurrence_total.l[,  "year" := years[j]]
      
      #what type of metric?
      dissimilarities_abundance_balanced.l[, "abund_biomass" := "biomass"]
      dissimilarities_occurrence_turnover.l[, "abund_biomass" := "biomass"]
      
      dissimilarities_abundance_gradient.l[, "abund_biomass" := "biomass"]
      dissimilarities_occurrence_nestedness.l[, "abund_biomass" := "biomass"]
      
      dissimilarities_abundance_total.l[, "abund_biomass" := "biomass"]
      dissimilarities_occurrence_total.l[, "abund_biomass" := "biomass"]
      
      #merge distance with dissimilarity for this year with both metrics of dissimilarity
      dissimilarities_full_wgt <- dissimilarities_abundance_balanced.l[dissimilarities_occurrence_turnover.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_wgt <- dissimilarities_full_wgt[dissimilarities_occurrence_total.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_wgt <- dissimilarities_full_wgt[dissimilarities_occurrence_nestedness.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_wgt <- dissimilarities_full_wgt[dissimilarities_abundance_gradient.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_wgt <- dissimilarities_full_wgt[dissimilarities_abundance_total.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      #add survey survey unit
      dissimilarities_full_wgt[, "survey" := scotland_shelf_sea.r_subset[1, survey]]
      dissimilarities_full_wgt[, "survey_unit" := scotland_survey_units[i]]
      
      #add to data table
      distances_dissimilarities_scotland_shelf_wgt  <- rbind(distances_dissimilarities_scotland_shelf_wgt ,  dissimilarities_full_wgt)
      
    ###ABUNDANCE
      
      
      #if some rows have num_cpue missing,  get rid of these rows
      reduced_year_abun <- reduced_year[complete.cases(reduced_year[, num_cpue_latband]), ]
      
      reduced_year_wide <- dcast(reduced_year_abun,  stat_rec_lat + year ~ accepted_name,  value.var = "num_cpue_latband",  fun.aggregate = sum) #longitude to wide data for community matrix,  column names are cell then species
      
      
      ncols <- ncol(reduced_year_wide)
      communitymatrix <- reduced_year_wide[, 3:ncols] #community matrix
      communitymatrix.occurence <- communitymatrix
      communitymatrix.occurence[communitymatrix.occurence > 0] <- 1
      
      dissimilarities_abundance <- beta.pair.abund(communitymatrix,  index.family = "bray") #dissimilarity 
      dissimilarities_occurrence <- beta.pair(communitymatrix.occurence,  index.family = "jaccard") #dissimilarity
      
      #make into matrix
      dissimilarities_abundance_balanced.m <- as.matrix(dissimilarities_abundance$beta.bray.bal,  labels=TRUE) #bal = balanced
      dissimilarities_abundance_gradient.m <- as.matrix(dissimilarities_abundance$beta.bray.gra,  labels=TRUE) #gra = gradient
      dissimilarities_abundance_total.m <- as.matrix(dissimilarities_abundance$beta.bray,  labels=TRUE) #total
      
      dissimilarities_occurrence_turnover.m <- as.matrix(dissimilarities_occurrence$beta.jtu,  labels=TRUE) #jtu = turnover
      dissimilarities_occurrence_nestedness.m <- as.matrix(dissimilarities_occurrence$beta.jne,  labels=TRUE) #jne = nestedness
      dissimilarities_occurrence_total.m <- as.matrix(dissimilarities_occurrence$beta.jac,  labels=TRUE) #total
      
      colnames(dissimilarities_abundance_balanced.m) <- rownames(dissimilarities_abundance_balanced.m) <- key_IDs_subset
      colnames(dissimilarities_abundance_gradient.m) <- rownames(dissimilarities_abundance_gradient.m) <- key_IDs_subset
      colnames(dissimilarities_abundance_total.m) <- rownames(dissimilarities_abundance_total.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_turnover.m) <- rownames(dissimilarities_occurrence_turnover.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_nestedness.m) <- rownames(dissimilarities_occurrence_nestedness.m) <- key_IDs_subset
      colnames(dissimilarities_occurrence_total.m) <- rownames(dissimilarities_occurrence_total.m) <- key_IDs_subset
      
      #reshape dissimilarities
      dissimilarities_abundance_balanced.l <- reshape2::melt(dissimilarities_abundance_balanced.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_balanced")
      dissimilarities_abundance_gradient.l <- reshape2::melt(dissimilarities_abundance_gradient.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_gradient")
      dissimilarities_abundance_total.l <- reshape2::melt(dissimilarities_abundance_total.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_total")
      
      dissimilarities_occurrence_turnover.l <- reshape2::melt(dissimilarities_occurrence_turnover.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_turnover")
      dissimilarities_occurrence_nestedness.l <- reshape2::melt(dissimilarities_occurrence_nestedness.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_nestedness")
      dissimilarities_occurrence_total.l <- reshape2::melt(dissimilarities_occurrence_total.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_total")
      
      #and then to data table format
      dissimilarities_abundance_balanced.l <- data.table(dissimilarities_abundance_balanced.l)
      dissimilarities_abundance_gradient.l <- data.table(dissimilarities_abundance_gradient.l)
      dissimilarities_abundance_total.l <- data.table(dissimilarities_abundance_total.l)
      dissimilarities_occurrence_turnover.l <- data.table(dissimilarities_occurrence_turnover.l)
      dissimilarities_occurrence_nestedness.l <- data.table(dissimilarities_occurrence_nestedness.l)
      dissimilarities_occurrence_total.l <- data.table(dissimilarities_occurrence_total.l)
      
      #add year for these values
      dissimilarities_abundance_balanced.l[,  "year" := years[j]]
      dissimilarities_occurrence_turnover.l[,  "year" := years[j]]
      
      
      dissimilarities_abundance_gradient.l[,  "year" := years[j]]
      dissimilarities_occurrence_nestedness.l[,  "year" := years[j]]
      
      dissimilarities_abundance_total.l[,  "year" := years[j]]
      dissimilarities_occurrence_total.l[,  "year" := years[j]]
      
      #what type of metric?
      dissimilarities_abundance_balanced.l[, "abund_biomass" := "abund"]
      dissimilarities_occurrence_turnover.l[, "abund_biomass" := "abund"]
      
      dissimilarities_abundance_gradient.l[, "abund_biomass" := "abund"]
      dissimilarities_occurrence_nestedness.l[, "abund_biomass" := "abund"]
      
      dissimilarities_abundance_total.l[, "abund_biomass" := "abund"]
      dissimilarities_occurrence_total.l[, "abund_biomass" := "abund"]
      
      #merge distance with dissimilarity for this year with both metrics of dissimilarity
      dissimilarities_full_abun <- dissimilarities_abundance_balanced.l[dissimilarities_occurrence_turnover.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_occurrence_total.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_occurrence_nestedness.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_abundance_gradient.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_abundance_total.l,  on = c("stat_rec_1",  "stat_rec_2",  "year",  "abund_biomass")]
      
      #add survey survey unit
      dissimilarities_full_abun[, "survey" := scotland_shelf_sea.r_subset[1, survey]]
      dissimilarities_full_abun[, "survey_unit" := scotland_survey_units[i]]
      
      #add to data table
      distances_dissimilarities_scotland_shelf_abun <- rbind(distances_dissimilarities_scotland_shelf_abun,  dissimilarities_full_abun)
      
  } #closes year
  
} #closes quarter

#merge two
distances_dissimilarities_scotland_shelf <- rbind(distances_dissimilarities_scotland_shelf_wgt, distances_dissimilarities_scotland_shelf_abun)

#delete if stat_rec_1 = stat_rec_2
nrow(distances_dissimilarities_scotland_shelf)
distances_dissimilarities_scotland_shelf <- distances_dissimilarities_scotland_shelf[stat_rec_1 != stat_rec_2,] #drops. #1400 1-1 comparisons
nrow(distances_dissimilarities_scotland_shelf)

#drop if we include both 1-2 and 2-1
distances_dissimilarities_scotland_shelf <- distances_dissimilarities_scotland_shelf[stat_rec_1 > stat_rec_2,] #drops duplicates, halves
nrow(distances_dissimilarities_scotland_shelf)

saveRDS(distances_dissimilarities_scotland_shelf,  here::here("output","dissimilarities", "distances_dissimilarities_scotland_shelf.rds"))

distances_dissimilarities_scotland_shelf <- readRDS(here::here("output","dissimilarities", "distances_dissimilarities_scotland_shelf.rds"))

#Flip to similarities


distances_dissimilarities_scotland_shelf[, jaccard_similarity_total := 1-jaccard_dissimilarity_total][,
                                        jaccard_similarity_turnover := 1-jaccard_dissimilarity_turnover][,
                                        jaccard_similarity_nestedness := 1-jaccard_dissimilarity_nestedness][,
                                        bray_curtis_similarity_total := 1-bray_curtis_dissimilarity_total][,
                                        bray_curtis_similarity_balanced := 1-bray_curtis_dissimilarity_balanced][,
                                        bray_curtis_similarity_gradient := 1-bray_curtis_dissimilarity_gradient]                                                               

#Average by year
#for similarities
distances_dissimilarities_scotland_shelf[, jaccard_similarity_total_avg := mean(jaccard_similarity_total),.(survey_unit, year,  abund_biomass)][,
                                         jaccard_similarity_turnover_avg := mean(jaccard_similarity_turnover),.(survey_unit, year,  abund_biomass)][,
                                          jaccard_similarity_nestedness_avg := mean(jaccard_similarity_nestedness),.(survey_unit, year,  abund_biomass)][,
                                       bray_curtis_similarity_total_avg := mean(bray_curtis_similarity_total),.(survey_unit, year,  abund_biomass)][,
                                        bray_curtis_similarity_balanced_avg := mean(bray_curtis_similarity_balanced),.(survey_unit, year,  abund_biomass)][,
                                       bray_curtis_similarity_gradient_avg := mean(bray_curtis_similarity_gradient),.(survey_unit, year,  abund_biomass)]

#for dissimilarities
distances_dissimilarities_scotland_shelf[, jaccard_dissimilarity_total_avg := mean(jaccard_dissimilarity_total),.(survey_unit, year,  abund_biomass)][,
                                                                                                                                                jaccard_dissimilarity_turnover_avg := mean(jaccard_dissimilarity_turnover),.(survey_unit, year,  abund_biomass)][,
                                                                                                                                                                                                                                                           jaccard_dissimilarity_nestedness_avg := mean(jaccard_dissimilarity_nestedness),.(survey_unit, year,  abund_biomass)][,
                                                                                                                                                                                                                                                                                                                                                                          bray_curtis_dissimilarity_total_avg := mean(bray_curtis_dissimilarity_total),.(survey_unit, year,  abund_biomass)][,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       bray_curtis_dissimilarity_balanced_avg := mean(bray_curtis_dissimilarity_balanced),.(survey_unit, year,  abund_biomass)][,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          bray_curtis_dissimilarity_gradient_avg := mean(bray_curtis_dissimilarity_gradient),.(survey_unit, year,  abund_biomass)]

#medians

distances_dissimilarities_scotland_shelf[, jaccard_similarity_total_median_allyears := median(jaccard_similarity_total),.(survey_unit,abund_biomass)][,
                            jaccard_similarity_turnover_median_allyears := median(jaccard_similarity_turnover),.(survey_unit,abund_biomass)][,
                             jaccard_similarity_nestedness_median_allyears := median(jaccard_similarity_nestedness),.(survey_unit,abund_biomass)][,
                            bray_curtis_similarity_total_median_allyears := median(bray_curtis_similarity_total),.(survey_unit,abund_biomass)][,
                             bray_curtis_similarity_balanced_median_allyears := median(bray_curtis_similarity_balanced),.(survey_unit,abund_biomass)][,
                            bray_curtis_similarity_gradient_median_allyears := median(bray_curtis_similarity_gradient),.(survey_unit,abund_biomass)]

#averages
distances_dissimilarities_scotland_shelf.unique <- unique(distances_dissimilarities_scotland_shelf[,.(jaccard_similarity_total_avg, jaccard_similarity_nestedness_avg, jaccard_similarity_turnover_avg,
                                                                                                      bray_curtis_similarity_balanced_avg, bray_curtis_similarity_gradient_avg, bray_curtis_similarity_total_avg,
                                                                                                      jaccard_dissimilarity_total_avg, jaccard_dissimilarity_nestedness_avg, jaccard_dissimilarity_turnover_avg,
                                                                                                      bray_curtis_dissimilarity_balanced_avg, bray_curtis_dissimilarity_gradient_avg, bray_curtis_dissimilarity_total_avg,
                                                                                                      survey_unit, year, abund_biomass)])


#Make plots to compare with Magurran
#Figure 2a (Jaccard similarity vs. year, boxplots)

ggplot(distances_dissimilarities_scotland_shelf[survey_unit == "SWC-IBTS-1"]) +
  geom_boxplot(aes(x = as.numeric(year),y=jaccard_similarity_total, group = year)) +
  facet_wrap(~survey_unit, ncol = 2) +
  lims(x = c(1986,2013)
       , y= c(0.4,0.75)
       ) +
  geom_smooth(aes(x = as.numeric(year),y=jaccard_similarity_total), method = "lm") +
 # geom_hline(yintercept = unique(distances_dissimilarities_scotland_shelf[survey_unit == "SWC-IBTS-1"]$jaccard_similarity_total_median_allyears), color= "salmon",  linewidth = 2) +
 labs(x= "Year", y = "Jaccard similarity") +
   theme_bw()

#2002-2004,2010,2013

ggplot(distances_dissimilarities_scotland_shelf[survey_unit == "SWC-IBTS-1"]) +
  geom_boxplot(aes(x = as.numeric(year),y=bray_curtis_similarity_total, group = year)) +
  facet_wrap(~survey_unit+abund_biomass, ncol = 1) +
  lims(x = c(1986,2013)
       , y= c(0.15,0.6)
       ) +
  geom_smooth(aes(x = as.numeric(year),y=bray_curtis_similarity_total), method = "lm") +
  labs(x= "Year", y = "Bray Curtis similarity") +
  theme_classic()

#Turnover and nestedness on one plot (showing dissimilarity instead of similarity) (like supplementary fig 5a)
ggplot(distances_dissimilarities_scotland_shelf.unique[survey_unit == "SWC-IBTS-1"]) +
  geom_point(aes(x = as.numeric(year),y=jaccard_dissimilarity_turnover_avg), color = "blue") +
  geom_point(aes(x = as.numeric(year),y=jaccard_dissimilarity_nestedness_avg), color = "forestgreen") +
  lims(x = c(1986,2013)) +
  labs(x= "Year", y = "Jaccard dissimilarity") +
  theme_bw()


  

########################################################################
#Now, let's try doing this but with data available online for Scotland  (Code partly from Juliano and Aurore in Fishglob repository)
########################################################################

#From Magurran et. al. 2015: Species and temperature data were downloaded from the ICES portal 
  #(DATRAS Fish Survey Data ‘Scottish West Coast Survey For Commercial Fish Species 1985–2013’ 
  #(Available at https://datras.ices.dk, accessed on 2014))."

#Pulled in from https://datras.ices.dk/Data_products/Download/Download_Data_public.aspx on August 22, 2023
#one file from early to 2010: "SWC-IBTS (Scottish West Coast Bottom Trawl Survey up to 2010)",
#another from 2010 forward: "SCOWCGFS (Scottish West Coast Groundfish Survey (from 2011))"


#not immediately sure what columns etc. will be

#Pre 2010:"/data/Scotland_magurran_comparison/scotland_raw_data_earlier_datras.csv"

#Post 2010: "/data/Scotland_magurran_comparison/scotland_raw_data_later_datras.csv"

#read in raw scotland data using DATRAS package
#metadata
hh.swc <- getDATRAS(record='HH', survey='SWC-IBTS', years=c(1985:2010), 
                    quarters=c(1:4)) #pre 2011
hh.scowcgfs <- getDATRAS(record='HH', survey='SCOWCGFS', years=c(2011:2020), 
                         quarters=c(1:4)) #post 2010

hh <- rbind(hh.swc,hh.scowcgfs)


#CPUE/length

hl.swc <- getDATRAS(record='HL', survey='SWC-IBTS', years=c(1985:2010), 
                    quarters=c(1:4)) #pre 2011
hl.scowcgfs <- getDATRAS(record='HL', survey='SCOWCGFS', years=c(2011:2020), 
                         quarters=c(1:4)) #post 2010

hl <- rbind(hl.swc,hl.scowcgfs)


##########################################################################################
#### CREATE A UNIQUE HAUL ID
##########################################################################################
hl$HaulID <- paste(hl$Survey, hl$Year,hl$Quarter, hl$Country, hl$Ship, hl$Gear, hl$StNo, 
                   hl$HaulNo)
hh$HaulID <- paste(hh$Survey, hh$Year,hh$Quarter, hh$Country, hh$Ship, hh$Gear, hh$StNo, 
                   hh$HaulNo)

# Is the HaulID unique?
hhn <- unique(hh$HaulID)
length(hhn)==nrow(hh)

# Only keep hauls where there is the length composition. 
hh <- subset(hh, hh$HaulID %in% hl$HaulID)
hl <- subset(hl, hl$HaulID %in% hh$HaulID)


##########################################################################################
#### MERGE HH and HL FILES
##########################################################################################

haulidhl <- sort(unique(hl$HaulID))
haulidhh <- sort(unique(hh$HaulID))
identical(haulidhh, haulidhl)
rm(haulidhh, haulidhl)

# remove some columns in hl
hl$SweepLngt <- hl$SpecCodeType <- hl$SpecCode <- hl$Sex <- hl$DateofCalculation <- NULL
hl$RecordType <- hl$GearEx <- NULL

# remove some columns in hh
hh$DateofCalculation <- hh$ThClineDepth <- hh$ThermoCline <- hh$SwellHeight <- NULL
hh$SwellDir <- hh$WindSpeed <- hh$WindDir <- hh$BotCurSpeed <- NULL
hh$BotCurDir <- hh$SurCurSpeed <- hh$SurCurDir <- hh$SpeedWater <- hh$TowDir <- NULL
hh$WgtGroundRope <- hh$KiteDim <- hh$Buoyancy <- hh$Tickler <- NULL
hh$DoorWgt <- hh$DoorSurface <- hh$WarpDen <- hh$Warpdia <- hh$Warplngt <- NULL
hh$Rigging  <- hh$HydroStNo <- hh$HaulLat <-  hh$HaulLong <- hh$DayNight <- NULL
hh$Stratum <- hh$TimeShot <- hh$Day <- hh$RecordType <- hh$GearExp <- hh$DoorType <- NULL

survey <- merge(hh, hl, by=c('HaulID','Survey','Quarter','Country','Ship','Gear','StNo','HaulNo','Year'), all.x=FALSE, all.y=TRUE)

nrow(survey)==nrow(hl)

survey <- survey %>% 
  dplyr::rename(SBT = BotTemp,
                SST = SurTemp,
                Speed = GroundSpeed,
                AphiaID = Valid_Aphia)

#check years included (should be 1985-2012 atleast)
unique(sort(survey$Year)) #correct

##########################################################################################
#### RESCALE DATA INTO ABUNDANCE FOR THE HAUL DURATION AND ABUNDANCE AT LENGTH (SKIP FOR NOW, JUST DO JACCARD)
##########################################################################################
# If Data Type=='C', abundance at length already readjusted with time so get back the 
# abundance for the actual duration of the haul.
# If data type=='R', abundance at length is multiplied by subfactor and adjusted to time
survey$CatCatchWgt = as.numeric(survey$CatCatchWgt)

#survey <- survey %>% 
#  mutate(HLNoAtLngt = case_when(DataType=='C' ~ HLNoAtLngt*SubFactor*HaulDur/60,
#                                DataType %in% c('S','R') ~ HLNoAtLngt*SubFactor),
#         TotalNo = case_when(DataType=='C' ~ TotalNo*HaulDur/60, 
#                             DataType %in% c('S','R') ~ TotalNo),
#         CatCatchWgt = case_when(DataType=='C' ~ CatCatchWgt*HaulDur/60,
#                                 DataType %in% c('S','R') ~ CatCatchWgt)) %>% 
#  select(-HaulVal, -DataType, -StdSpecRecCode, -SpecVal, -SubWgt, -SubFactor) %>% 
#  mutate(Survey = if_else(Survey=='SCOWCGFS', 'SWC-IBTS', Survey)) %>% 
#  mutate(Survey = if_else(Survey=='SCOROC','ROCKALL',Survey)) %>% 
#  filter(!(Survey=="NS-IBTS" & BySpecRecCode %in% c(0,2,3,4,5)), 
#         # remove hauls where not all species are recorded
#         !(Survey=="BITS" & BySpecRecCode==0))


#Magurran et al. 2015 steps
  # Species abundance = CPUE (# individuals per species per 1h tow)
  # 35 rectangles with good community time series (colname = StatRec)
ICES_rec_include <- c("47E3","47E4","47E5",
                      "46E2","46E3","46E4","46E5",
                      "45E0","45E1","45E2","45E3","45E4",
                      "44E0","44E1","44E3","44E4",
                      "43E1","43E3",
                      "42E0","42E1","42E2","42E3",
                      "41E0","41E1","41E2","41E3",
                      "40E0","40E1","40E2","40E3",
                      "39E0","39E1","39E2","39E4","39E5")
      #vs
  unique(survey$StatRec)
  #-->
    survey.dt <- data.table(survey)
    survey.dt <- survey.dt[StatRec %in% ICES_rec_include,]
  #131 taxa sampled, plus 5 macroinverts, total of 126 finfish
  length(unique(survey.dt$AphiaID)) #this dataset includes 178, but we'll keep for now
  
  #1992 and 1995 have substantially fewer samples, so they are excluded
  #-->
  survey.dt <- survey.dt[!(Year %in% c(1992,1995)),]
    
  #Trawling occurs in the first quarter of the year (January, February and March)
    #-->
    survey.dt <- survey.dt[Month %in% c(1,2,3),]
  
  #Assign rectangles to nine 30' latitudinal bands
    survey.dt[,stat_rec_lat := substr(StatRec,1,2)]
    #TotalNo = Total number of fish in the given haul and species
    survey.dt[,num_cpue_latband :=
                         sum(TotalNo, na.rm = T),.(stat_rec_lat,Year, AphiaID)]
    
  #compile community time series for each latitudinal band
    survey.dt.r <- unique(survey.dt[,.(Year, stat_rec_lat, AphiaID, num_cpue_latband)])
    
    fwrite(survey.dt.r, file.path("data","scotland_magurran_comparison","survey.dt.r.csv")) #now, 160 unique species
    
  #they conduct sample rarefaction, but we will skip this step for now as it's not relevant for spatial beta diversity because we must maintain species IDs
    
###########################################
### Perform occurrence and abundance based spatial beta diversity dissimilarity calculations over time with this 'raw' dataset straight from DATRAS
###########################################

    #loop through all years
    
distances_dissimilarities_scotland_shelf_DATRAS <- data.table("stat_rec_1" = integer(), 
                                                                "stat_rec_2" = integer(), 
                                                                "bray_curtis_dissimilarity_balanced" = numeric(), 
                                                                year = integer(), 
                                                                "jaccard_dissimilarity_turnover" = numeric(), 
                                                                "jaccard_dissimilarity_nestedness" =  numeric(), 
                                                                "bray_curtis_dissimilarity_gradient" = as.numeric(), 
                                                                "jaccard_dissimilarity_total" = numeric(), 
                                                                "bray_curtis_dissimilarity_total" = numeric())
      
years <- sort(unique(survey.dt.r[, as.numeric(Year)]))
      
      #haul id keys
lat_recs <- survey.dt.r[,.(stat_rec_lat)]
      
    for (j in 1:length(years)) {
        reduced_year <- survey.dt.r[Year == years[j], ]
        
        ###ABUNDANCE (not going through length weight conversion to do this for weight)
        
        
        #if some rows have num_cpue missing,  get rid of these rows
        reduced_year_abun <- reduced_year[complete.cases(reduced_year[, num_cpue_latband]), ]
        
        reduced_year_wide <- dcast(reduced_year_abun,  stat_rec_lat + Year ~ AphiaID,  value.var = "num_cpue_latband",  fun.aggregate = sum) #longitude to wide data for community matrix,  column names are cell then species
        
        #list of stat_rec_lats
        key_IDs_subset <- reduced_year_wide[,stat_rec_lat]
        
        ncols <- ncol(reduced_year_wide)
        communitymatrix <- reduced_year_wide[, 3:ncols] #community matrix
        communitymatrix.occurence <- communitymatrix
        communitymatrix.occurence[communitymatrix.occurence > 0] <- 1
        
        dissimilarities_abundance <- beta.pair.abund(communitymatrix,  index.family = "bray") #dissimilarity 
        dissimilarities_occurrence <- beta.pair(communitymatrix.occurence,  index.family = "jaccard") #dissimilarity
        
        #make into matrix
        dissimilarities_abundance_balanced.m <- as.matrix(dissimilarities_abundance$beta.bray.bal,  labels=TRUE) #bal = balanced
        dissimilarities_abundance_gradient.m <- as.matrix(dissimilarities_abundance$beta.bray.gra,  labels=TRUE) #gra = gradient
        dissimilarities_abundance_total.m <- as.matrix(dissimilarities_abundance$beta.bray,  labels=TRUE) #total
        
        dissimilarities_occurrence_turnover.m <- as.matrix(dissimilarities_occurrence$beta.jtu,  labels=TRUE) #jtu = turnover
        dissimilarities_occurrence_nestedness.m <- as.matrix(dissimilarities_occurrence$beta.jne,  labels=TRUE) #jne = nestedness
        dissimilarities_occurrence_total.m <- as.matrix(dissimilarities_occurrence$beta.jac,  labels=TRUE) #total
        
        colnames(dissimilarities_abundance_balanced.m) <- rownames(dissimilarities_abundance_balanced.m) <- key_IDs_subset
        colnames(dissimilarities_abundance_gradient.m) <- rownames(dissimilarities_abundance_gradient.m) <- key_IDs_subset
        colnames(dissimilarities_abundance_total.m) <- rownames(dissimilarities_abundance_total.m) <- key_IDs_subset
        colnames(dissimilarities_occurrence_turnover.m) <- rownames(dissimilarities_occurrence_turnover.m) <- key_IDs_subset
        colnames(dissimilarities_occurrence_nestedness.m) <- rownames(dissimilarities_occurrence_nestedness.m) <- key_IDs_subset
        colnames(dissimilarities_occurrence_total.m) <- rownames(dissimilarities_occurrence_total.m) <- key_IDs_subset
        
        #reshape dissimilarities
        dissimilarities_abundance_balanced.l <- reshape2::melt(dissimilarities_abundance_balanced.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_balanced")
        dissimilarities_abundance_gradient.l <- reshape2::melt(dissimilarities_abundance_gradient.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_gradient")
        dissimilarities_abundance_total.l <- reshape2::melt(dissimilarities_abundance_total.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "bray_curtis_dissimilarity_total")
        
        dissimilarities_occurrence_turnover.l <- reshape2::melt(dissimilarities_occurrence_turnover.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_turnover")
        dissimilarities_occurrence_nestedness.l <- reshape2::melt(dissimilarities_occurrence_nestedness.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_nestedness")
        dissimilarities_occurrence_total.l <- reshape2::melt(dissimilarities_occurrence_total.m,  varnames = c("stat_rec_1",  "stat_rec_2"),  value.name = "jaccard_dissimilarity_total")
        
        #and then to data table format
        dissimilarities_abundance_balanced.l <- data.table(dissimilarities_abundance_balanced.l)
        dissimilarities_abundance_gradient.l <- data.table(dissimilarities_abundance_gradient.l)
        dissimilarities_abundance_total.l <- data.table(dissimilarities_abundance_total.l)
        dissimilarities_occurrence_turnover.l <- data.table(dissimilarities_occurrence_turnover.l)
        dissimilarities_occurrence_nestedness.l <- data.table(dissimilarities_occurrence_nestedness.l)
        dissimilarities_occurrence_total.l <- data.table(dissimilarities_occurrence_total.l)
        
        #add year for these values
        dissimilarities_abundance_balanced.l[,  "year" := years[j]]
        dissimilarities_occurrence_turnover.l[,  "year" := years[j]]
        
        
        dissimilarities_abundance_gradient.l[,  "year" := years[j]]
        dissimilarities_occurrence_nestedness.l[,  "year" := years[j]]
        
        dissimilarities_abundance_total.l[,  "year" := years[j]]
        dissimilarities_occurrence_total.l[,  "year" := years[j]]
        
        #merge distance with dissimilarity for this year with both metrics of dissimilarity
        dissimilarities_full_abun <- dissimilarities_abundance_balanced.l[dissimilarities_occurrence_turnover.l,  on = c("stat_rec_1",  "stat_rec_2",  "year")]
        
        dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_occurrence_total.l,  on = c("stat_rec_1",  "stat_rec_2",  "year")]
        
        dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_occurrence_nestedness.l,  on = c("stat_rec_1",  "stat_rec_2",  "year")]
        
        dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_abundance_gradient.l,  on = c("stat_rec_1",  "stat_rec_2",  "year")]
        
        dissimilarities_full_abun <- dissimilarities_full_abun[dissimilarities_abundance_total.l,  on = c("stat_rec_1",  "stat_rec_2",  "year")]
        
        #add to data table
        distances_dissimilarities_scotland_shelf_DATRAS <- rbind(distances_dissimilarities_scotland_shelf_DATRAS,  dissimilarities_full_abun)
        
      } #closes year

    
    #delete if stat_rec_1 = stat_rec_2
    nrow(distances_dissimilarities_scotland_shelf_DATRAS)
    distances_dissimilarities_scotland_shelf_DATRAS <- distances_dissimilarities_scotland_shelf_DATRAS[stat_rec_1 != stat_rec_2,] #drops 1-1 comparisons
    nrow(distances_dissimilarities_scotland_shelf_DATRAS)
    
    #drop if we include both 1-2 and 2-1
    distances_dissimilarities_scotland_shelf_DATRAS <- distances_dissimilarities_scotland_shelf_DATRAS[stat_rec_1 > stat_rec_2,] #drops duplicates, halves
    nrow(distances_dissimilarities_scotland_shelf_DATRAS) #1224
    
    saveRDS(distances_dissimilarities_scotland_shelf_DATRAS,  here::here("output","dissimilarities", "distances_dissimilarities_scotland_shelf_DATRAS.rds"))
    
    distances_dissimilarities_scotland_shelf_DATRAS <- readRDS(here::here("output","dissimilarities", "distances_dissimilarities_scotland_shelf_DATRAS.rds"))
    
    #Flip to similarities
    
    
    distances_dissimilarities_scotland_shelf_DATRAS[, jaccard_similarity_total := 1-jaccard_dissimilarity_total][,
                                                                                                          jaccard_similarity_turnover := 1-jaccard_dissimilarity_turnover][,
                                                                                                                                                                           jaccard_similarity_nestedness := 1-jaccard_dissimilarity_nestedness][,
                                                                                                                                                                                                                                                bray_curtis_similarity_total := 1-bray_curtis_dissimilarity_total][,
                                                                                                                                                                                                                                                                                                                   bray_curtis_similarity_balanced := 1-bray_curtis_dissimilarity_balanced][,
                                                                                                                                                                                                                                                                                                                                                                                            bray_curtis_similarity_gradient := 1-bray_curtis_dissimilarity_gradient]                                                               
    
    #Average by year
    
    #similarities
    distances_dissimilarities_scotland_shelf_DATRAS[, jaccard_similarity_total_avg := mean(jaccard_similarity_total),.(year)][,
                                                                                                                                                    jaccard_similarity_turnover_avg := mean(jaccard_similarity_turnover),.(year)][,
                                                                                                                                                                                                                                                               jaccard_similarity_nestedness_avg := mean(jaccard_similarity_nestedness),.(year)][,
                                                                                                                                                                                                                                                                                                                                                                              bray_curtis_similarity_total_avg := mean(bray_curtis_similarity_total),.(year)][,bray_curtis_similarity_gradient_avg := mean(bray_curtis_similarity_gradient),.(year)]
 #by dissimilarities
    distances_dissimilarities_scotland_shelf_DATRAS[, jaccard_dissimilarity_total_avg := mean(jaccard_dissimilarity_total),.(year)][,
                                                                                                                              jaccard_dissimilarity_turnover_avg := mean(jaccard_dissimilarity_turnover),.(year)][,
                                                                                                                                                                                                            jaccard_dissimilarity_nestedness_avg := mean(jaccard_dissimilarity_nestedness),.(year)][,
                                                                                                                                                                                                                                                                                              bray_curtis_dissimilarity_total_avg := mean(bray_curtis_dissimilarity_total),.(year)][,
                                                                                                                                                                                                                                                                                                                                                                              bray_curtis_dissimilarity_balanced_avg := mean(bray_curtis_dissimilarity_balanced),.(year)][,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    bray_curtis_dissimilarity_gradient_avg := mean(bray_curtis_dissimilarity_gradient),.(year)]
    
    
    
    
    #averages
    distances_dissimilarities_scotland_shelf_DATRAS.unique <- unique(distances_dissimilarities_scotland_shelf_DATRAS[,.(year, jaccard_similarity_total_avg, jaccard_similarity_nestedness_avg, jaccard_similarity_turnover_avg,
                                                                                                          bray_curtis_similarity_balanced_avg, bray_curtis_similarity_gradient_avg, bray_curtis_similarity_total_avg,
                                                                                                          jaccard_dissimilarity_total_avg, jaccard_dissimilarity_nestedness_avg, jaccard_dissimilarity_turnover_avg,
                                                                                                          bray_curtis_dissimilarity_balanced_avg, bray_curtis_dissimilarity_gradient_avg, bray_curtis_dissimilarity_total_avg)])
    
    
    #Again, make plots to compare with Magurran
    #Figure 2a (Jaccard similarity vs. year, boxplots)
    
    ggplot(distances_dissimilarities_scotland_shelf_DATRAS) +
      geom_boxplot(aes(x = as.numeric(year),y=jaccard_similarity_total, group = year)) +
      lims(x = c(1986,2013)
           , y= c(0.4,0.75)
      ) +
      geom_smooth(aes(x = as.numeric(year),y=jaccard_similarity_total), method = "lm") +
      # geom_hline(yintercept = unique(distances_dissimilarities_scotland_shelf[survey_unit == "SWC-IBTS-1"]$jaccard_similarity_total_median_allyears), color= "salmon",  linewidth = 2) +
      labs(x= "Year", y = "Jaccard similarity") +
      theme_bw()
    
#Turnover and nestedness on one plot (showing dissimilarity instead of similarity) (like supplementary fig 5a)
    ggplot(distances_dissimilarities_scotland_shelf_DATRAS.unique) +
      geom_point(aes(x = as.numeric(year),y=jaccard_dissimilarity_turnover_avg), color = "blue") +
      geom_point(aes(x = as.numeric(year),y=jaccard_dissimilarity_nestedness_avg), color = "forestgreen") +
      lims(x = c(1986,2013)) +
      labs(x= "Year", y = "Jaccard dissimilarity") +
      theme_bw()
    
#TAKEAWAY
    #While DATRAS doesn't match our FISHGLOB dataset perfectly, the patterns are more similar than the patterns of Magurran et al. 2015. Therefore, I will move forward with our dataset.
    