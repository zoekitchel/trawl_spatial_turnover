#Greenland exploration


#Working with Karl Michael Werner to try to figure out why the lowest 15% species doesn't align with his expectations for 
#Greenland bottom trawl survey data. From this examination, we found out that Greenland data have errors, and therefore 
#pull in separate Greenland data outside of Fishglob v1.5.


library(data.table)
library(ggplot2)
library(sf)
library(rgeos)
library(concaveman)
library(raster)
library(stringr)
library(rnaturalearth)
library(rmapshaper)
library(nlme)
library(MuMIn)
library(gridExtra)


#raw data
FishGlob_clean <- readRDS(here::here("data", "cleaned", "FishGlob.wellsampledyearscells_complete.final.rds"))

#delete all rows where wgt_cpue is NA,  unless it is MEDITS region
FishGlob_clean.noNA <- FishGlob_clean[(survey != "MEDITS" & !is.na(wgt_cpue)) | survey == "MEDITS"]

#check
summary(FishGlob_clean.noNA[, wgt_cpue]) #why inf? 13 observations from GSL-N,  I will delete these
FishGlob_clean.noNA <- copy(FishGlob_clean.noNA)[!is.infinite(wgt_cpue), ]

FishGlob_clean.noNA[, wgt := sum(wgt,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                               wgt_h := sum(wgt_h,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                            wgt_cpue := sum(wgt_cpue,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                               num := sum(num,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                                                                                        num_h := sum(num_h,  na.rm = T), .(accepted_name,  haul_id)][, 
                                                                                                                                                                                                                                                                                                                                     num_cpue := sum(num_cpue,  na.rm = T), .(accepted_name,  haul_id)]

FishGlob_clean.noNA <- unique(FishGlob_clean.noNA[,.(survey, haul_id, country, sub_area, continent, stat_rec, station, stratum, year, month, day, quarter, season, latitude, longitude, haul_dur, area_swept, gear, depth, sbt, sst, num, num_h, num_cpue, wgt, wgt_h, wgt_cpue, verbatim_aphia_id, accepted_name, aphia_id, SpecCode, kingdom, phylum, class, order, family, genus, rank, survey_unit, years_sampled, cell_ID, cell_year_count, year_cell_count, years_sampled_update)])

#make depth numeric
FishGlob_clean.noNA[,depth := as.numeric(depth)]

#make survey and survey unit factors
FishGlob_clean.noNA[,survey:=factor(survey)][,survey_unit:=factor(survey_unit)]

#adjust years
FishGlob_clean.noNA[,year_adj := year-min(year)+1]

#only greenland
Greenland <- FishGlob_clean.noNA[survey_unit == "GRL-DE"]

#sum wgt for each year
Greenland[,wgt_cpue_sum := sum(wgt_cpue, na.rm = T),.(year,accepted_name)]

#Greenland summed by year
Greenland.r <- unique(Greenland[,.(accepted_name, year, wgt_cpue_sum)])

#Independent time series for all species
# Create a list to store plots
plot_list <- list()

#names of species
GRL_spp <- levels(as.factor(Greenland.r$accepted_name))

# Generate some sample data and create plots
for (i in 1:length(GRL_spp)) {
  #timeseries for single species
  Greenland_onespp <- Greenland.r[accepted_name == GRL_spp[i]]
  
  # Create a ggplot object
  p <- ggplot(Greenland_onespp, aes(year, wgt_cpue_sum)) +
    geom_point() +
    labs(x = "Year",y="CPUE in kg/km^2") +
    ggtitle(paste(GRL_spp[i])) +
    ylim(0,NA) +
    xlim(1993,2015) +
      theme_classic()
  
  plot_list[[i]] <- p
}


# Create a function to arrange plots in a single column and save to a PDF document
save_plots_to_single_column_pdf <- function(plots, filename) {
  num_plots <- length(plots)
  pdf(file = filename, height = 11, width = 8.5)  # Set PDF dimensions (adjust as needed)
  
  for (page_start in seq(from = 1, to = num_plots, by = 5)) {
    page_end <- min(page_start + 4, num_plots)
    grid.arrange(grobs = plots[page_start:page_end], ncol = 1)
  }
  
  dev.off()
}

# Call the function to save plots to a single-column PDF document
save_plots_to_single_column_pdf(plot_list, "figures/Greenland_spp_plots.pdf")

