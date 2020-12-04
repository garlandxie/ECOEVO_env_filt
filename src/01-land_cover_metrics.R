################################################################################
# Accompanying code for the paper: 
#   No environmental filtering of wasps and bees during urbanization
#
# Paper authors: 
#   Garland Xie          (1), 
#   Nicholas Sookhan     (1), 
#   Kelly Carscadden     (2), 
#   and J. Scott MacIvor (1)
#
# Corresponding authors for this script:  
#   Garland Xie      (1)
#   Nicholas Sookhan (2)
#
# Affiliations: 
#   (1) Department of Biological Sciences, 
#       University of Toronto Scarborough,
#       1265 Military Trail, Toronto, ON, M1C 1A4, Canada
#       email: garland.xie@mail.utoronto.ca, 
#              nicholas.sookhan@mail.utoronto.ca
#              scott.macivor@mail.utoronto.ca
#   (2) Department of Ecology and Evolutionary Biology,
#       University of Colorado Boulder,
#
# Purpose of this R script: to clean up data for the land cover metrics 
#
# IMPORTANT: Running this code is a bit slow (approx. 15 minutes).

# libraries --------------------------------------------------------------------
library(here)              # for creating relative file-paths
library(raster)            # for manipulating raster datasets
library(sf)                # for manipulating GIS data
library(readxl)            # for importing excel files
library(landscapemetrics)  # for calculating landscape composition metrics
library(dplyr)             # for manipulating data frames
library(tidyr)             # for pivoting tables from long to wide format
library(rgdal)
# import -----------------------------------------------------------------------

# land cover data
lc <- raster(
  here("data/original", 
       "toronto_2007_landcover.img")
)

# abundance per site 
site<- read.csv(
  here("data/original", "site.csv"), 
  row.names = 1
)

# functions
source(here("src", "functions.R"))

# data clean -------------------------------------------------------------------

# get landcover raster  projection
lc_proj <- proj4string(lc)

# ensure that site data is in the same coordinate system as the raster data
point <- st_as_sf(site, coords = c("Longitude", "Latitude"), crs = 4326)
point <- st_transform(point,lc_proj)

# buffers: 250 spatial scale ---------------------------------------------------

# create 250m buffer radii
buffer_250 <- st_buffer(point, 250)

# convert from a data frame to a list 
buffer_250 <- split(buffer_250[,1], f = buffer_250[,1,drop = T])

#only retain buffers which intersect lc 
overl_250 <- unlist(
  lapply(
    buffer_250, 
    function(x) 
      class(raster::intersect(
        extent(lc),
        x)) == "Extent"
   )  
  )

# get insersecting buffers
buffer_250 <- buffer_250[overl_250]

# buffers: 500m spatial scale --------------------------------------------------

# create 500m buffer radii
buffer_500 <- st_buffer(point, 500)

# convert from a data frame to a list 
buffer_500 <- split(buffer_500[,1], f = buffer_500[,1,drop=T])

# only retain buffers which intersect lc 
overl_500 <- unlist(lapply(buffer_500, 
                           function(x) class(raster::intersect(extent(lc),x))=="Extent"
                           )
                   )
# get intersecting buffers
buffer_500 <- buffer_500[overl_500]

# landscape composition: 250 spatial scale -------------------------------------

pland_250 <- list(rep(NA, times = length(buffer_250)))

for (k in 1:length(buffer_250)) {
  pland_250[[k]] <- calc_pland(l = lc, buffer = buffer_250[[k]])
}

land_use_250 <- do.call("rbind", pland_250) 

# landscape composition: 500m spatial scale ------------------------------------

pland_500 <- list(rep(NA, times = length(buffer_500)))

for (i in 1:length(buffer_500)) {
    pland_500[[i]] <- calc_pland(l = lc, buffer = buffer_500[[i]])
}

land_use_500 <- do.call("rbind", pland_500)

# calc missing data: 250 -------------------------------------------------------

prop_miss_250 <- list(rep(NA, times = length(buffer_250)))

for (k in 1:length(buffer_250)) {
  prop_miss_250[[k]] <- calc_prop_miss(l = lc, buffer = buffer_250[[k]])
}

prop_miss_250 <- do.call("rbind", prop_miss_250)

# sites

# 102 = "nSPob" (boundary edge) water
# 36 = "BW4HD" (boundary edge) water
# 166 - "XvuFi" (odd) water
# 58 - fLwK1 (boundary)
# 87 - LOMxN (boundary edge)
# 156 -  vtkfe (boundary edge)
# 47 - E3dsm (boundary edge)
# 90 - Lwe2a (boundary edge)
# 75 - iMm0w (boundary)
# 38 - CAgNl (boundary)
# 41 - Ci7F8 (boundary)
# 73 - IBsVs (boundary)
# 105 - OUdi4 (odd)

# calc missing data: 500 -------------------------------------------------------

prop_miss_500 <- list(rep(NA, times = length(buffer_500)))

for (k in 1:length(buffer_500)) {
  prop_miss_500[[k]] <- calc_prop_miss(l = lc, buffer = buffer_500[[k]])
}

prop_miss_500 <- do.call("rbind", prop_miss_500)
  
# clean: land cover ------------------------------------------------------------

# 250m
lw_250 <- land_use_250 %>% 
  select(class, id, prop_land_use) %>%
  pivot_wider(names_from = class, values_from = prop_land_use) %>% 
  select(site = `id`,
         prop_tree_250   = `1`,
         prop_grass_250  = `2`,
         prop_earth_250  = `3`,
         prop_water_250  = `4`,
         prop_build_250  = `5`,
         prop_roads_250  = `6`,
         prop_paved_250  = `7`,
         prop_agri_250   = `8`
  ) %>%
  mutate(prop_urb_250 = prop_roads_250 + prop_paved_250 + prop_build_250,
         across(where(is.numeric), ~replace_na(., 0)))

# 500m
lw_500 <- land_use_500 %>% 
  select(class, id, prop_land_use) %>%
  pivot_wider(names_from = class, values_from = prop_land_use) %>% 
  select(site = `id`,
         prop_tree_500 = `1`,
         prop_grass_500 = `2`,
         prop_earth_500 = `3`,
         prop_water_500 = `4`,
         prop_build_500 = `5`,
         prop_roads_500 = `6`,
         prop_paved_500 = `7`,
         prop_agri_500  = `8`
  ) %>%
  mutate(prop_urb_500 = prop_roads_500 + prop_paved_500 + prop_build_500,
         across(where(is.numeric), ~replace_na(., 0)))

# clean: remove sites  ---------------------------------------------------------
outside_TO <- c(
  "GAJVv", 
  "SQWq3",
  "N53op", 
  "auCMf", 
  "lWpWV",
  "Z42dv"
)

# remove sites that are outside raster boundaries 
l_250 <- lw_250 %>%
    filter(!site %in% outside_TO)

l_500 <- lw_500 %>%
    filter(!site %in% outside_TO)

# save to disk -----------------------------------------------------------------
write.csv(l_250, 
          file = here(
            "data/working",
            "land_use_250.csv")
          )

write.csv(l_500,
          file = here(
            "data/working", 
            "land_use_500.csv")
          )
