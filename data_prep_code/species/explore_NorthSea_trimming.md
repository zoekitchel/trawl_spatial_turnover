Exploring North Sea data trimming
================

  - [Load data](#load-data)
  - [Standardize trawl footprint](#standardize-trawl-footprint)
      - [Hex size 8](#hex-size-8)
          - [Plot num hauls per cell and
            year](#plot-num-hauls-per-cell-and-year)
          - [Trim out years with \<70% cells in a
            survey](#trim-out-years-with-70-cells-in-a-survey)
          - [Plot num years by cell](#plot-num-years-by-cell)
          - [Mark cells sampled \>=3x each remaining
            year](#mark-cells-sampled-3x-each-remaining-year)
      - [Hex size 7](#hex-size-7)
          - [Plot num hauls per cell and
            year](#plot-num-hauls-per-cell-and-year-1)
          - [Trim out years with \<70% cells in a
            survey](#trim-out-years-with-70-cells-in-a-survey-1)
          - [Plot num years by cell](#plot-num-years-by-cell-1)
          - [Mark cells sampled \>=3x each remaining
            year](#mark-cells-sampled-3x-each-remaining-year-1)

``` r
library(dggridR)
```

    ## Loading required package: dplyr

    ## Warning: package 'dplyr' was built under R version 4.0.5

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: rlang

    ## Warning: package 'rlang' was built under R version 4.0.5

    ## Loading required package: sf

    ## Warning: package 'sf' was built under R version 4.0.5

    ## Linking to GEOS 3.9.1, GDAL 3.4.0, PROJ 8.1.1; sf_use_s2() is TRUE

    ## Loading required package: sp

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:rlang':
    ## 
    ##     :=

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
library(rgdal)
```

    ## rgdal: version: 1.5-23, (SVN revision 1121)
    ## Geospatial Data Abstraction Library extensions to R successfully loaded
    ## Loaded GDAL runtime: GDAL 3.2.1, released 2020/12/29
    ## Path to GDAL shared files: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rgdal/gdal
    ## GDAL binary built with GEOS: TRUE 
    ## Loaded PROJ runtime: Rel. 7.2.1, January 1st, 2021, [PJ_VERSION: 721]
    ## Path to PROJ shared files: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rgdal/proj
    ## PROJ CDN enabled: FALSE
    ## Linking to sp version:1.4-5
    ## To mute warnings of possible GDAL/OSR exportToProj4() degradation,
    ## use options("rgdal_show_exportToProj4_warnings"="none") before loading rgdal.
    ## Overwritten PROJ_LIB was /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rgdal/proj

``` r
library(raster)
```

    ## 
    ## Attaching package: 'raster'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     shift

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(sp)
library(here)
```

    ## here() starts at /Users/mpinsky/Documents/Rutgers/Community_and_climate/trawl_spatial_turnover

``` r
library(sf)
library(ggplot2) # for plotting
library(gridExtra) # for arranging plots
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

# Load data

Pull in compiled and cleaned data from fishglob February 4, 2022  
Trim to North Sea and divide into quarters 1 and 3

``` r
FishGlob <- fread(file = here::here("temp", "FISHGLOB_v1.1_clean.csv"))
FishGlob <- FishGlob[survey == "NS-IBTS",]
FishGlob[,survey_season := as.factor(paste0(survey,"_",quarter))] # quarter indicates season
FishGlob[, summary(survey_season)]
```

    ## NS-IBTS_1 NS-IBTS_3 
    ##    228975    135189

``` r
#if positive, subtract 360
FishGlob[,longitude_s := ifelse(longitude > 150,(longitude-360),(longitude))]

#delete if NA for longitude or latitude
FishGlob <- FishGlob[complete.cases(FishGlob[,.(longitude, latitude)])]
```

# Standardize trawl footprint

## Hex size 8

``` r
#if positive, subtract 360
FishGlob[,longitude_s := ifelse(longitude > 150,(longitude-360),(longitude))]

#delete if NA for longitude or latitude
FishGlob <- FishGlob[complete.cases(FishGlob[,.(longitude, latitude)])]

#set up grid
dggs <- dgconstruct(res = 8, metric = T) #with res = 8, we will need at least 3 observations per year within 7,774.2 km^2 (roughly size of some NEUS strata)

#pull out unique lat lons and get grid cells
unique_latlon <- unique(FishGlob[,.(latitude, longitude_s)])
unique_latlon[,cell := dgGEO_to_SEQNUM(dggs, longitude_s, latitude)] #get corresponding grid cells for this region/survey combo

#find cell centers
cellcenters <- dgSEQNUM_to_GEO(dggs, unique_latlon[,cell])

#linking cell centers to unique_latlon
unique_latlon[,cell_center_longitude_s := cellcenters$lon_deg][,cell_center_latitude:= cellcenters$lat_deg]

#link centers back to main data table
FishGlob.dg <- merge(FishGlob, unique_latlon, by = c("latitude", "longitude_s"), all.x = TRUE)
```

### Plot num hauls per cell and year

Red lines maximize the number of cell x years without any missing.  
Purple lines maximize the number of cell x years with \<5% missing.

``` r
# make a list of all unique year x cell x survey
year_cell_count.dt <- as.data.table(expand.grid(cell = unique(FishGlob.dg$cell), 
                                                year = unique(FishGlob.dg$year), 
                                                survey_season = unique(FishGlob.dg$survey_season))) 

# calc num hauls per year x cell x survey, then merge to full list
temp <- FishGlob.dg[, .(nhaul = length(unique(haul_id))), by = .(cell, year, survey_season)]
year_cell_count.dt <- merge(year_cell_count.dt, temp, all.x = TRUE)
year_cell_count.dt[is.na(nhaul), nhaul := 0] # fill in for year x cell x survey that aren't present

# loop through each survey x season. Couldn't figure out how to do this without a loop.
survseas <- FishGlob.dg[, unique(survey_season)]
plots <- vector('list', length(survseas))
for(i in 1:length(survseas)){
  # trim to this survey
  thisyearcellcount <- year_cell_count.dt[survey_season == survseas[i],]
  
  # calculate the most to least sampled cells and years
  year_order <- thisyearcellcount[nhaul >0, .(ncell = length(unique(cell))), 
                                     by = year][ncell >0,][order(ncell),]
  cell_order <- thisyearcellcount[nhaul > 0, .(nyear = length(unique(year))), 
                                     by = cell][nyear >0,][order(nyear),]
  
  # set the factor order by ave number of hauls
  thisyearcellcount[, year := factor(year, levels = year_order$year)]
  thisyearcellcount[, cell := factor(cell, levels = cell_order$cell)]
  
  # calculate num missing cells x years for different thresholds. slow.
  cutoffs <- data.table(expand.grid(year = levels(thisyearcellcount$year), cell = levels(thisyearcellcount$cell)))
  cutoffs[, ':='(ntot = NA_integer_, nmiss = NA_integer_, nkeep = NA_integer_)]
  nyr <- length(unique(cutoffs$year))
  ncl <- length(unique(cutoffs$cell))
  for(j in 1:nyr){
    for(k in 1:ncl){
      thisntot <- thisyearcellcount[cell %in% levels(cell)[k:ncl] & year %in% levels(year)[j:nyr], .N]
      thisnmiss <- thisyearcellcount[cell %in% levels(cell)[k:ncl] & year %in% levels(year)[j:nyr], sum(nhaul == 0)]
      thisnkeep <- thisyearcellcount[cell %in% levels(cell)[k:ncl] & year %in% levels(year)[j:nyr], sum(nhaul > 0)]
      cutoffs[cell == thisyearcellcount[, levels(cell)[k]] &
                year == thisyearcellcount[,levels(year)[j]],
              ':='(ntot = thisntot, nmiss = thisnmiss, nkeep = thisnkeep)]
    }
  }
  
  # choose a threshold 
  chosencutoffs0 <- cutoffs[nmiss==0,][nkeep == max(nkeep),] # based on nothing missing
  chosencutoffs02 <- cutoffs[nmiss/ntot < 0.02,][nkeep == max(nkeep),] # based on <2% missing
  
  # make plot
  plots[[i]] <- ggplot(thisyearcellcount[nhaul > 0,], aes(year, cell, color = nhaul)) +
    geom_point() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
          axis.text.y = element_text(size = 5)) +
    labs(title = survseas[i]) +
    geom_vline(xintercept = chosencutoffs0[, as.numeric(year)], color = 'red') +
    geom_hline(yintercept = chosencutoffs0[, as.numeric(cell)], color = 'red') +
    geom_vline(xintercept = chosencutoffs02[, as.numeric(year)], color = 'purple') +
    geom_hline(yintercept = chosencutoffs02[, as.numeric(cell)], color = 'purple')
}

# plots
do.call("grid.arrange", c(plots, ncol=2))
```

![](explore_NorthSea_trimming_files/figure-gfm/plot%20num%20hauls%20per%20cell%20and%20year%20hex%208-1.png)<!-- -->

### Trim out years with \<70% cells in a survey

yearkeep70 column == FALSE for years to trim out

``` r
# Find max num cells in each survey, calc 70% treshold
year_cells <- FishGlob.dg[,.(cell_count = length(unique(cell))), by = .(year, survey_season)]
year_cells[,thresh := 0.7*(max(cell_count)), by = survey_season] # of cells/ year to cut off below
  
#only keep years where over 70% of cells are sampled
year_cells[,yearkeep70 := cell_count > thresh]
  
# mark years to keep in FishGlob (yearkeep70 colum)
FishGlob.yr <- merge(FishGlob.dg, year_cells[, .(year, survey_season, yearkeep70)], all.x = TRUE, by = c('year', 'survey_season'))
```

#### Plot num. hauls by cell latitude and year.

Second plot shows whether year is kept or not.

``` r
# table of years by cells, count is # hauls
yearcell_table <- FishGlob.yr[, .(nhaul = length(unique(haul_id))), by = .(cell, year, cell_center_latitude, cell_center_longitude_s, yearkeep70, survey_season)]

# color by #hauls
ggplot(yearcell_table, aes(year, cell_center_latitude, color = nhaul)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# color by keep or not at 70% threshold
ggplot(yearcell_table, aes(year, cell_center_latitude, color = yearkeep70)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

### Plot num years by cell

Only count years from those remaining after the previous trim Second
plot shows whether a cell is sampled in all years

``` r
# table of years by cells, count is # hauls
cell_table <- FishGlob.yr[yearkeep70 == TRUE, .(nyear = length(unique(year))), by = .(cell, cell_center_latitude, cell_center_longitude_s, survey_season)]
cell_table[, allyrs := nyear == max(nyear), by = survey_season]

ggplot(cell_table, aes(cell_center_longitude_s, cell_center_latitude, color = nyear)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(cell_table, aes(cell_center_longitude_s, cell_center_latitude, color = allyrs)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

### Mark cells sampled \>=3x each remaining year

``` r
#identify num hauls per cells x year
FishGlob.yr[,year_cell_count := length(unique(haul_id)),.(year,cell, survey_season)]

#num years in which cells are sampled at least X times in any remaining year
# use X=1 as a test: should match plot of cells sampled in all years
nyrs_by_cell1dt <- FishGlob.yr[year_cell_count >= 1 & yearkeep70, .(nyrs_by_cell1 = length(unique(year))), by = .(cell, survey_season)] # number of years in which this cell has >= 1 haul
nyrs_by_cell3dt <- FishGlob.yr[year_cell_count >= 3 & yearkeep70, .(nyrs_by_cell3 = length(unique(year))), by = .(cell, survey_season)] # number of years in which this cell has >= 1 haul

FishGlob.cell <- merge(FishGlob.yr, nyrs_by_cell1dt, all.x = TRUE, by = c('survey_season', 'cell'))
FishGlob.cell[is.na(nyrs_by_cell1)==TRUE, nyrs_by_cell1 := 0] # fill missing values
FishGlob.cell <- merge(FishGlob.cell, nyrs_by_cell3dt, all.x = TRUE, by = c('survey_season', 'cell'))
FishGlob.cell[is.na(nyrs_by_cell3), nyrs_by_cell3 := 0] # fill missing values

# mark cells sampled the max years
FishGlob.cell[, cellkeep1 := (nyrs_by_cell1 == length(unique(year))) & yearkeep70, by = .(survey_season, yearkeep70)] # TRUE if all remaining years met the min num hauls threshold (1 in this case)
FishGlob.cell[, cellkeep3 := (nyrs_by_cell3 == length(unique(year))) & yearkeep70, by = .(survey_season, yearkeep70)]

#make a map of these points
ggplot(FishGlob.cell[yearkeep70 == TRUE][!duplicated(cbind(cell, survey_season)), ], 
       aes(cell_center_longitude_s, cell_center_latitude, color = cellkeep3)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Hex size 7

``` r
#set up grid
dggs7 <- dgconstruct(res = 7, metric = T) # larger than res 8

#pull out unique lat lons and get grid cells
unique_latlon7 <- unique(FishGlob[,.(latitude, longitude_s)])
unique_latlon7[,cell := dgGEO_to_SEQNUM(dggs7, longitude_s, latitude)] #get corresponding grid cells for this region/survey combo

#find cell centers
cellcenters7 <- dgSEQNUM_to_GEO(dggs7, unique_latlon7[,cell])

#linking cell centers to unique_latlon
unique_latlon7[,cell_center_longitude_s := cellcenters7$lon_deg][,cell_center_latitude:= cellcenters7$lat_deg]

#link centers back to main data table
FishGlob.dg7 <- merge(FishGlob, unique_latlon7, by = c("latitude", "longitude_s"), all.x = TRUE)
```

### Plot num hauls per cell and year

Red lines maximize the number of cell x years without any missing.  
Purple lines maximize the number of cell x years with \<5% missing.

``` r
# make a list of all unique year x cell x survey
year_cell_count.dt7 <- as.data.table(expand.grid(cell = unique(FishGlob.dg7$cell), 
                                                year = unique(FishGlob.dg7$year), 
                                                survey_season = unique(FishGlob.dg7$survey_season))) 

# calc num hauls per year x cell x survey, then merge to full list
temp7 <- FishGlob.dg7[, .(nhaul = length(unique(haul_id))), by = .(cell, year, survey_season)]
year_cell_count.dt7 <- merge(year_cell_count.dt7, temp7, all.x = TRUE)
year_cell_count.dt7[is.na(nhaul), nhaul := 0] # fill in for year x cell x survey that aren't present

# loop through each survey x season
survseas <- FishGlob.dg7[, unique(survey_season)]
plots <- vector('list', length(survseas))
for(i in 1:length(survseas)){
  # trim to this survey
  thisyearcellcount <- year_cell_count.dt7[survey_season == survseas[i],]
  
  # calculate the most to least sampled cells and years
  year_order <- thisyearcellcount[nhaul >0, .(ncell = length(unique(cell))), 
                                     by = year][ncell >0,][order(ncell),]
  cell_order <- thisyearcellcount[nhaul > 0, .(nyear = length(unique(year))), 
                                     by = cell][nyear >0,][order(nyear),]
  
  # set the factor order by ave number of hauls
  thisyearcellcount[, year := factor(year, levels = year_order$year)]
  thisyearcellcount[, cell := factor(cell, levels = cell_order$cell)]
  
  # calculate num missing cells x years for different thresholds. slow.
  cutoffs <- data.table(expand.grid(year = levels(thisyearcellcount$year), cell = levels(thisyearcellcount$cell)))
  cutoffs[, ':='(ntot = NA_integer_, nmiss = NA_integer_, nkeep = NA_integer_)]
  nyr <- length(unique(cutoffs$year))
  ncl <- length(unique(cutoffs$cell))
  for(j in 1:nyr){
    for(k in 1:ncl){
      thisntot <- thisyearcellcount[cell %in% levels(cell)[k:ncl] & year %in% levels(year)[j:nyr], .N]
      thisnmiss <- thisyearcellcount[cell %in% levels(cell)[k:ncl] & year %in% levels(year)[j:nyr], sum(nhaul == 0)]
      thisnkeep <- thisyearcellcount[cell %in% levels(cell)[k:ncl] & year %in% levels(year)[j:nyr], sum(nhaul > 0)]
      cutoffs[cell == thisyearcellcount[, levels(cell)[k]] &
                year == thisyearcellcount[,levels(year)[j]],
              ':='(ntot = thisntot, nmiss = thisnmiss, nkeep = thisnkeep)]
    }
  }
  
  # choose a threshold 
  chosencutoffs0 <- cutoffs[nmiss==0,][nkeep == max(nkeep),] # based on nothing missing
  chosencutoffs02 <- cutoffs[nmiss/ntot < 0.02,][nkeep == max(nkeep),] # based on <2% missing
  
  # make plot
  plots[[i]] <- ggplot(thisyearcellcount[nhaul > 0,], aes(year, cell, color = nhaul)) +
    geom_point() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
          axis.text.y = element_text(size = 5)) +
    labs(title = survseas[i]) +
    geom_vline(xintercept = chosencutoffs0[, as.numeric(year)], color = 'red') +
    geom_hline(yintercept = chosencutoffs0[, as.numeric(cell)], color = 'red') +
    geom_vline(xintercept = chosencutoffs02[, as.numeric(year)], color = 'purple') +
    geom_hline(yintercept = chosencutoffs02[, as.numeric(cell)], color = 'purple')
}

# plots
do.call("grid.arrange", c(plots, ncol=2))
```

![](explore_NorthSea_trimming_files/figure-gfm/plot%20num%20hauls%20per%20cell%20and%20year%20hex%207-1.png)<!-- -->

### Trim out years with \<70% cells in a survey

yearkeep70 column == FALSE for years to trim out

``` r
# Find max num cells in each survey, calc 70% treshold
year_cells7 <- FishGlob.dg7[,.(cell_count = length(unique(cell))), by = .(year, survey_season)]
year_cells7[,thresh := 0.7*(max(cell_count)), by = survey_season] # of cells/ year to cut off below
  
#only keep years where over 70% of cells are sampled
year_cells7[,yearkeep70 := cell_count > thresh]
  
# mark years to keep in FishGlob (yearkeep70 colum)
FishGlob.yr7 <- merge(FishGlob.dg7, year_cells7[, .(year, survey_season, yearkeep70)], all.x = TRUE, by = c('year', 'survey_season'))
```

#### Plot num. hauls by cell latitude and year.

Second plot shows whether year is kept or not.

``` r
# table of years by cells, count is # hauls
yearcell_table7 <- FishGlob.yr7[, .(nhaul = length(unique(haul_id))), by = .(cell, year, cell_center_latitude, cell_center_longitude_s, yearkeep70, survey_season)]

# color by #hauls
ggplot(yearcell_table7, aes(year, cell_center_latitude, color = nhaul)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# color by keep or not at 70% threshold
ggplot(yearcell_table7, aes(year, cell_center_latitude, color = yearkeep70)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

### Plot num years by cell

Only count years from those remaining after the previous trim Second
plot shows whether a cell is sampled in all years

``` r
# table of years by cells, count is # hauls
cell_table7 <- FishGlob.yr7[yearkeep70 == TRUE, .(nyear = length(unique(year))), by = .(cell, cell_center_latitude, cell_center_longitude_s, survey_season)]
cell_table7[, allyrs := nyear == max(nyear), by = survey_season]

ggplot(cell_table7, aes(cell_center_longitude_s, cell_center_latitude, color = nyear)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggplot(cell_table7, aes(cell_center_longitude_s, cell_center_latitude, color = allyrs)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

### Mark cells sampled \>=3x each remaining year

``` r
#identify num hauls per cells x year
FishGlob.yr7[,year_cell_count := length(unique(haul_id)),.(year,cell, survey_season)]

#num years in which cells are sampled at least X times in any remaining year
# use X=1 as a test: should match plot of cells sampled in all years
nyrs_by_cell1dt7 <- FishGlob.yr7[year_cell_count >= 1 & yearkeep70, .(nyrs_by_cell1 = length(unique(year))), by = .(cell, survey_season)] # number of years in which this cell has >= 1 haul
nyrs_by_cell3dt7 <- FishGlob.yr7[year_cell_count >= 3 & yearkeep70, .(nyrs_by_cell3 = length(unique(year))), by = .(cell, survey_season)] # number of years in which this cell has >= 1 haul

FishGlob.cell7 <- merge(FishGlob.yr7, nyrs_by_cell1dt7, all.x = TRUE, by = c('survey_season', 'cell'))
FishGlob.cell7[is.na(nyrs_by_cell1)==TRUE, nyrs_by_cell1 := 0] # fill missing values
FishGlob.cell7 <- merge(FishGlob.cell7, nyrs_by_cell3dt7, all.x = TRUE, by = c('survey_season', 'cell'))
FishGlob.cell7[is.na(nyrs_by_cell3), nyrs_by_cell3 := 0] # fill missing values

# mark cells sampled the max years
FishGlob.cell7[, cellkeep1 := (nyrs_by_cell1 == length(unique(year))) & yearkeep70, by = .(survey_season, yearkeep70)] # TRUE if all remaining years met the min num hauls threshold (1 in this case)
FishGlob.cell7[, cellkeep3 := (nyrs_by_cell3 == length(unique(year))) & yearkeep70, by = .(survey_season, yearkeep70)]

#make a map of these points
ggplot(FishGlob.cell7[yearkeep70 == TRUE][!duplicated(cbind(cell, survey_season)), ], 
       aes(cell_center_longitude_s, cell_center_latitude, color = cellkeep3)) +
  geom_point() +
  facet_grid(~survey_season)
```

![](explore_NorthSea_trimming_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
