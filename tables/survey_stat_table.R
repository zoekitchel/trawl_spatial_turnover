library(data.table)
library(vegan)
library(sf)
library(concaveman) #polygon around points
library(geosphere)
library(ggpubr) #stat_regline_equation
library(nlme)


FishGlob_clean <- readRDS(here::here("data", "cleaned", "FishGlob.wellsampledyearscells_complete.final.rds"))


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



#Calculate survey details only

survey_stats <- data.table(survey = character() ,  
                           survey_unit = character(),  spp_num = numeric(),  
                           study_period = numeric(),  study_duration = numeric(),  
                           lat_range = numeric(),  mid_lat = numeric(),  lon_range = numeric(),  
                           area = numeric(),  depth_range = numeric(),   mid_depth = numeric())


distances_dissimilarities_allyears <- data.table("haul_id1" = integer(), 
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

for (i in 1:length(all_survey_units)) {
  
  FishGlob_clean_subset <- FishGlob_clean.noNA[survey_unit == all_survey_units[i], ]
  
  #map
  ####unique lat lon
  #add column with adjusted longitude for few surveys that cross dateline (NZ-CHAT and AI)
  FishGlob_clean_subset[,longitude_adj := ifelse((survey_unit %in% c("AI","NZ-CHAT") & longitude > 0),longitude-360,longitude)]
  
  lat_lon <- unique(FishGlob_clean_subset[, .(latitude,  longitude_adj)])
  
  pts <- st_as_sf(lat_lon,  coords=c('longitude_adj', 'latitude'),  crs=4326 )
  
  conc <- concaveman(pts,  concavity = 5)
  sf::sf_use_s2(FALSE) #turn off the s2 processing; GEOS treats projected coordinates as planar (i.e. two points lie on a line of infinite max lenght) while s2 is more "correct" (two points lie on a great circle of circumference of 40 075 kilometers)
  area <- st_area(conc) #m2,  check this later
  
  #fill row
  row <- data.table(FishGlob_clean_subset[1, survey], 
                    FishGlob_clean_subset[1, as.character(survey_unit)], 
                    as.numeric(length(unique(FishGlob_clean_subset[, accepted_name]))), 
                    as.numeric(max(FishGlob_clean_subset[, year])-min(FishGlob_clean_subset[, year])), 
                    as.numeric(length(unique(FishGlob_clean_subset[, year]))), 
                    as.numeric(max(FishGlob_clean_subset[, latitude])-min(FishGlob_clean_subset[, latitude])), 
                    as.numeric(mean(FishGlob_clean_subset[, latitude])), 
                    as.numeric(max(FishGlob_clean_subset[, longitude_adj])-min(FishGlob_clean_subset[, longitude_adj])), 
                    as.numeric(area), 
                    as.numeric(max(FishGlob_clean_subset[, as.numeric(depth)],  na.rm = T)-min(FishGlob_clean_subset[, as.numeric(depth)],  na.rm = T)), 
                    as.numeric(mean(FishGlob_clean_subset[, as.numeric(depth)],  na.rm = T)))
  
  survey_stats <- rbind(survey_stats,  row,  use.names = F)
  
  
}

saveRDS(survey_stats, file.path("tables", "survey_stats.Rds"))

#TO DO
#clean survey stats
