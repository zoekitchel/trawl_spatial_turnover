#only dissimilariities to speed up process

#SAME AS annotate_dissim_script_3.r, but only for 85% most abundant species

library(data.table)
library(vegan)
library(betapart) #allows us to partition beta diversity
library(geosphere)
library(nlme)
library(reshape2)



FishGlob_clean <- readRDS(here::here("homogenization_analyses", "FishGlob.wellsampledyearscells_complete.final.rds"))



#delete all rows where wgt_cpue is NA,  unless it is MEDITS region
FishGlob_clean.noNA <- FishGlob_clean[(survey != "MEDITS" & !is.na(wgt_cpue)) | survey == "MEDITS"]

#check
summary(FishGlob_clean.noNA[, wgt_cpue]) #why inf? 13 observations from GSL-N,  I will delete these
FishGlob_clean.noNA <- copy(FishGlob_clean.noNA)[!is.infinite(wgt_cpue), ]

all_survey_units <- unique(FishGlob_clean.noNA[, survey_unit])


FishGlob_clean.noNA[survey_unit == "NZ-CHAT" & month == 12,  year := year+1, ]



FishGlob_clean.noNA[, wgt := sum(wgt,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                               wgt_h := sum(wgt_h,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                            wgt_cpue := sum(wgt_cpue,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                               num := sum(num,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                                                                                        num_h := sum(num_h,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                                                                                                                                                     num_cpue := sum(num_cpue,  na.rm = T), .(accepted_name,  haul_id)]

FishGlob_clean.noNA <- unique(FishGlob_clean.noNA[,.(survey, haul_id, country, sub_area, continent, stat_rec, station, stratum, year, month, day, quarter, season, latitude, longitude, haul_dur, area_swept, gear, depth, sbt, sst, num, num_h, num_cpue, wgt, wgt_h, wgt_cpue, verbatim_aphia_id, accepted_name, aphia_id, SpecCode, kingdom, phylum, class, order, family, genus, rank, survey_unit, years_sampled, cell_ID, cell_year_count, year_cell_count, years_sampled_update)])


#add new columns that sum biomass for each individual species per survey, and 
FishGlob_clean.noNA[survey_unit != "MEDITS",abundance_spp_sum := sum(wgt_cpue),.(accepted_name, survey_unit)][survey_unit != "MEDITS",abundance_sum := sum(wgt_cpue),survey_unit][survey_unit == "MEDITS",abundance_spp_sum := sum(num_cpue),.(accepted_name, survey_unit)][survey_unit == "MEDITS",abundance_sum := sum(num_cpue),survey_unit]



FishGlob_abundance_survey <- unique(FishGlob_clean.noNA[,.(survey_unit, accepted_name, abundance_spp_sum, abundance_sum)])

FishGlob_abundance_survey[,perc_spp_total_abundance := abundance_spp_sum/abundance_sum][,spp_count := uniqueN(accepted_name),survey_unit][,spp_count15_percent := round(0.15*spp_count,0)]

#rank 15% least abundant spp 
#rank by abundance
setkey(FishGlob_abundance_survey, perc_spp_total_abundance)

FishGlob_abundance_survey[, rank_abundance_bysurvey := frank(perc_spp_total_abundance), survey_unit]

#code each spp as 90% most common by abundance or 15% lease common by abundance in a survey
FishGlob_abundance_survey[, exclude_spp := ifelse(rank_abundance_bysurvey <= spp_count15_percent, T,F)]

#spp and survey key
spp_survey_keep_key <- unique(FishGlob_abundance_survey[,.(survey_unit, accepted_name, exclude_spp)])

FishGlob_clean.noNA_15perc <- spp_survey_keep_key[FishGlob_clean.noNA, on = .(survey_unit, accepted_name)]

FishGlob_clean.noNA_15perc_excluded <- FishGlob_clean.noNA_15perc_omit[exclude_spp == F,]

distances_dissimilarities_allyears_15perc_excluded <- data.table("haul_id1" = integer(), 
                                                 "haul_id2" = integer(), 
                                                 "distance" = numeric(), 
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

#Dissimilarity loop
for (i in 1:length(all_survey_units)) {
  FishGlob_clean_subset <- FishGlob_clean.noNA_15perc_excluded[survey_unit == all_survey_units[i], ]
  
  #map
  ####unique lat lon
  #add column with adjusted longitude for few surveys that cross dateline (NZ-CHAT and AI)
  FishGlob_clean_subset[,longitude_adj := ifelse((survey_unit %in% c("AI","NZ-CHAT") & longitude > 0),longitude-360,longitude)]
  
  lat_lon <- unique(FishGlob_clean_subset[, .(latitude,  longitude_adj)])
  
  #list years
  FishGlob_clean_subset[, year:= as.numeric(year)] #make numeric
  setorder(FishGlob_clean_subset,  year)
  years <- unique(FishGlob_clean_subset[, year])
  
  #haul id keys
  haul_ids <- unique(FishGlob_clean_subset[, haul_id])
  haul_ids_key <- data.table(haul_id = haul_ids,  key_ID = seq(1, length(haul_ids),  by = 1))
  
  
  #convert haul_ids to numeric key_ids
  FishGlob_clean_subset <- FishGlob_clean_subset[haul_ids_key,  on = "haul_id"]
  
  for (j in 1:length(years)) {
    reduced_year <- FishGlob_clean_subset[year == years[j], ]
    
    #distances among cells
    setorder(reduced_year,  key_ID)
    
    latitude_longitude_haul_id <- unique(reduced_year[, .(latitude, longitude, key_ID)])
    distances <- distm(latitude_longitude_haul_id[, .(longitude, latitude)])
    key_IDs <- latitude_longitude_haul_id[, key_ID]
    
    colnames(distances) <- rownames(distances) <- key_IDs
    
    #wide to long
    haul_id_distances.l <- reshape2::melt(distances, varnames = (c("haul_id1",  "haul_id2")),  value.name = "distance")
    
    #make into data table
    haul_id_distances.l <- data.table(haul_id_distances.l)
    
    if(!(reduced_year[1, survey] == "MEDITS")) {
      
      #if some rows have wgt_cpue missing,  get rid of these rows
      reduced_year <- reduced_year[complete.cases(reduced_year[, wgt_cpue]), ]
      
      reduced_year_wide <- dcast(reduced_year,  key_ID + year ~ accepted_name,  value.var = "wgt_cpue",  fun.aggregate = sum) #longitude to wide data for community matrix,  column names are cell then species
      
      
      ncols <- ncol(reduced_year_wide)
      communitymatrix <- reduced_year_wide[, 3:ncols] #community matrix
      communitymatrix.occurence <- communitymatrix
      communitymatrix.occurence[communitymatrix.occurence > 0] <- 1
      
      #list of haul_id keys
      key_IDs_subset <- reduced_year_wide$key_ID
      
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
      dissimilarities_abundance_balanced.l <- reshape2::melt(dissimilarities_abundance_balanced.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_balanced")
      dissimilarities_abundance_gradient.l <- reshape2::melt(dissimilarities_abundance_gradient.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_gradient")
      dissimilarities_abundance_total.l <- reshape2::melt(dissimilarities_abundance_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_total")
      
      dissimilarities_occurrence_turnover.l <- reshape2::melt(dissimilarities_occurrence_turnover.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_turnover")
      dissimilarities_occurrence_nestedness.l <- reshape2::melt(dissimilarities_occurrence_nestedness.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_nestedness")
      dissimilarities_occurrence_total.l <- reshape2::melt(dissimilarities_occurrence_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_total")
      
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
      dissimilarities_full <- haul_id_distances.l[dissimilarities_abundance_balanced.l,  on = c("haul_id1",  "haul_id2")]
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_turnover.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_total.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_nestedness.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_gradient.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_total.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      #add survey survey unit
      dissimilarities_full[, "survey" := FishGlob_clean_subset[1, survey]]
      dissimilarities_full[, "survey_unit" := all_survey_units[i]]
      
      #add to data table
      distances_dissimilarities_allyears_15perc_excluded <- rbind(distances_dissimilarities_allyears_15perc_excluded,  dissimilarities_full)
      
      print(paste0(j, "/", length(years)))
      
    } else { #if we do using abundance instead
      
      #if some rows have num_cpue missing,  get rid of these rows
      reduced_year <- reduced_year[complete.cases(reduced_year[, num_cpue]), ]
      
      reduced_year_wide <- dcast(reduced_year,  key_ID + year ~ accepted_name,  value.var = "num_cpue",  fun.aggregate = sum) #longitude to wide data for community matrix,  column names are cell then species
      
      
      ncols <- ncol(reduced_year_wide)
      communitymatrix <- reduced_year_wide[, 3:ncols] #community matrix
      communitymatrix.occurence <- communitymatrix
      communitymatrix.occurence[communitymatrix.occurence > 0] <- 1
      
      #list of haul_id keys
      key_IDs_subset <- reduced_year_wide$key_ID
      
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
      dissimilarities_abundance_balanced.l <- reshape2::melt(dissimilarities_abundance_balanced.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_balanced")
      dissimilarities_abundance_gradient.l <- reshape2::melt(dissimilarities_abundance_gradient.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_gradient")
      dissimilarities_abundance_total.l <- reshape2::melt(dissimilarities_abundance_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "bray_curtis_dissimilarity_total")
      
      dissimilarities_occurrence_turnover.l <- reshape2::melt(dissimilarities_occurrence_turnover.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_turnover")
      dissimilarities_occurrence_nestedness.l <- reshape2::melt(dissimilarities_occurrence_nestedness.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_nestedness")
      dissimilarities_occurrence_total.l <- reshape2::melt(dissimilarities_occurrence_total.m,  varnames = c("haul_id1",  "haul_id2"),  value.name = "jaccard_dissimilarity_total")
      
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
      dissimilarities_full <- haul_id_distances.l[dissimilarities_abundance_balanced.l,  on = c("haul_id1",  "haul_id2")]
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_turnover.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_total.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_occurrence_nestedness.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_gradient.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      dissimilarities_full <- dissimilarities_full[dissimilarities_abundance_total.l,  on = c("haul_id1",  "haul_id2",  "year",  "abund_biomass")]
      
      #add survey survey unit
      dissimilarities_full[, "survey" := FishGlob_clean_subset[1, survey]]
      dissimilarities_full[, "survey_unit" := all_survey_units[i]]
      
      #add to data table
      distances_dissimilarities_allyears_15perc_excluded <- rbind(distances_dissimilarities_allyears_15perc_excluded,  dissimilarities_full)
      
      print(paste0(j, "/", length(years)))
      
    }
    
    
    
    
  } 
  
  print(paste0("We have cycled through ",  i, " of ", length(all_survey_units),  " total survey/season combos"))
  
}

#closes survey/region


saveRDS(distances_dissimilarities_allyears_15perc_excluded,  here::here("homogenization_analyses", "distances_dissimilarities_allyears_15perc_excluded.rds"))
