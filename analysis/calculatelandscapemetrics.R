##### Get landscape metrics at site level
### April 14, 2026
### J Melanson


### Load in packages
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)
library(terra)
library(raster)
library(sf)
library(landscapemetrics)

bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"

#####################################################
# Load and floral survey data
#####################################################
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T))
veg2022 = read.csv(paste0(bombus_path, "raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "raw_data/2023vegetationdata.csv"))
sample2022 = read.csv(paste0(bombus_path, "raw_data/2022sampledata.csv"))
sample2023 = read.csv(paste0(bombus_path, "raw_data/2023sampledata.csv"))


#####################################################
# Load and prep landcover data
#####################################################
landscape_raster = raster(paste0(bombus_path, "landscape/full_rasters/FValley_2res_32610.tif"))
fv_points = st_read(paste0(bombus_path, "landscape/fvbombus/fvbombus_points.shp"))
fv_points$site_name[fv_points$site_name == "pit_meadows"] = "pitt_meadows"

#Change CRS to meters
fv_points = st_transform(fv_points, 32610)
landscape_raster = rast(landscape_raster)
crs(landscape_raster) = "EPSG:32610"


# reclassification matrix for iji
# seminatural == forest, wetland, grassland (5, 14, 17)
# out of production == blackberry, grassy margin, hedgerow, fallow (1, 4, 7, 9)
# annual crops = annual, polyculture (0, 13)
# perennial crops == blueberry, cranberry, other perennial (2, 3, 12)
# managed grassy == hay, pasture (6, 11)
# urban = urban (15)
# water = marine, water (10, 16)
# industrial = industrial (8)

rcl_matrix = rbind(
  cbind(c(0, 13), 1),        # annual crops
  cbind(c(1, 7, 9), 2),  # field margins
  cbind(4, 3),           # fallow
  cbind(5,4),   # forest
  cbind(14,5), # seminatural grassland
  cbind(17,6), # wetlands
  cbind(2, 7), # blueberry
  cbind(3, 8), # cranberry
  cbind(12, 9),    # other perennials
  cbind(c(6, 11), 10),        # managed grassy
  cbind(15, 11),              # urban
  cbind(c(10, 16), 12),      # water
  cbind(8, 13)                # industrial
)

compute_site_metrics <- function(site_points, raster, buffer_dist = 500) {
  
  # buffer all points in the site and dissolve into one polygon
  site_buffer = site_points %>%
    st_buffer(buffer_dist) %>%
    st_union()
  
  # crop and mask raster to the buffered area
  site_vect   = vect(site_buffer)
  site_raster = crop(raster, site_vect) %>% mask(site_vect)
  
  # reclassify and crop for iji
  site_raster_rcl = classify(site_raster, rcl = rcl_matrix) %>%
    mask(site_vect)
  
  # composition metrics on original raster
  comp_metrics = calculate_lsm(
    site_raster,
    what = c("lsm_c_pland"))
  
  # configuration metric on reclassified raster
  config_metrics = calculate_lsm(
    site_raster_rcl,
    what = "lsm_l_iji")
  
  # bind together and return
  bind_rows(comp_metrics, config_metrics)
}

# apply for all six sites
allmetrics = fv_points %>%
  group_by(site_name) %>%
  group_modify(~ compute_site_metrics(.x, landscape_raster)) %>%
  bind_rows()

iji = allmetrics %>% filter(metric == "iji")
iji = iji[c("site_name", "metric", "value")]
iji$year = NA

comp = allmetrics %>% filter(metric != "iji")

#######################################
# Compute flowering habitat prop
#######################################
# floral habitat classes = 
# blackberry (1) ##### SEASON DEPENDENT (2022: 157-233 , 2023: 146-234)
# blueberry (2) ##### SEASON DEPENDENT (2022: 130-174, 2023: 114-162)
# cranberry (3) ##### SEASON DEPENDENT (from https://www.bccranberries.com/wp-content/uploads/2022/03/BCCRF_2021_January2022_FinalReport-harbut.pdf -- June and July e.g., 152-212)
# fallow (4)
# forest (5)
# hedgerow (7)
# lowedge (9)
# seminatgrass (14)
# urban (15)
# wetland (17)

# compute 7 different groupings: 
# none, blueberry, blue+cran, blue+black, all, black+cran, black
# plus permanent habitats

# no berries
none = comp %>%
  filter(class %in% c(4, 5, 7, 9, 14, 15, 17)) %>%
  group_by(site_name) %>%
  summarize(value = sum(value))

# with blueberries
blue = comp %>%
  filter(class %in% c(4, 5, 7, 9, 14, 15, 17, 2)) %>%
  group_by(site_name) %>%
  summarize(value = sum(value))

# with blueberries and cranberries
bluecran = comp %>%
  filter(class %in% c(4, 5, 7, 9, 14, 15, 17, 2, 3)) %>%
  group_by(site_name) %>%
  summarize(value = sum(value))

# with blueberries and blackberries
blueblack = comp %>%
  filter(class %in% c(4, 5, 7, 9, 14, 15, 17, 2, 1)) %>%
  group_by(site_name) %>%
  summarize(value = sum(value))

# with all three
allberry = comp %>%
  filter(class %in% c(4, 5, 7, 9, 14, 15, 17, 1, 2, 3)) %>%
  group_by(site_name) %>%
  summarize(value = sum(value))

# with cranberry and blackberry
cranblack = comp %>%
  filter(class %in% c(4, 5, 7, 9, 14, 15, 17, 1, 3)) %>%
  group_by(site_name) %>%
  summarize(value = sum(value))

# with blackberry
black = comp %>%
  filter(class %in% c(4, 5, 7, 9, 14, 15, 17, 1)) %>%
  group_by(site_name) %>%
  summarize(value = sum(value))



#######################################
# Make temporally explicit
#######################################

floralpheno2022 = bind_rows(
  blue  %>% crossing(julian_date = 129:151),
  bluecran %>% crossing(julian_date = 152:156),
  allberry %>% crossing(julian_date = 157:173),
  cranblack %>% crossing(julian_date = 174:212),
  black  %>% crossing(julian_date = 213:233))
floralpheno2022$year = 2022
floralpheno2022$metric = "floweringpercent"

floralpheno2023 = bind_rows(
  none %>% crossing(julian_date = 81:113),
  blue  %>% crossing(julian_date = 114:145),
  blueblack %>% crossing(julian_date = 146:151),
  allberry %>% crossing(julian_date = 152:162),
  cranblack %>% crossing(julian_date = 163:212),
  black  %>% crossing(julian_date = 213:233))
floralpheno2023$year = 2023
floralpheno2023$metric = "floweringpercent"

alllandscapemetrics = rbind(iji, floralpheno2022, floralpheno2023)
write.csv(alllandscapemetrics, "analysis/landscapemetrics.csv")
