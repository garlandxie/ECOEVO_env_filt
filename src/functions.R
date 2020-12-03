# Calculate landscape composition metrics
#
# @description 
# This function calculates landscape composition 
# metrics from a specific buffer.

# * `crop()` clips the raster to the extent of a given buffer
# * `rasterize()` ensures that the clipped object has the geometry of a buffer
# * `lsm_c_pland()` calculates the landscape composition metrics
#
# @return a data frame of land cover class ID, site ID, and
# proportion of a given land cover class
#
# @param l A RasterLayer object that represents the 2008 
# Forest and Land Cover raster dataset from Open Data Toronto
# 
# @param buffer A list of sf objects of the buffers for all the sampled sites
#
# @details 
# Link to raster data: https://open.toronto.ca/dataset/forest-and-land-cover/
# Note that the landcover classes are coded as follows: 
# (1) tree canopy, 
# (2) grass/shrub, 
# (3) bare earth, 
# (4) water, 
# (5) buildings, 
# (6) roads, 
# (7) other paved surfaces, 
# (8) agriculture 

calc_Stats <- function(l, buffer) {
  
  clip1 <- raster::crop(l, extent(buffer))
  clip2 <- raster::rasterize(buffer, clip1, mask = TRUE)
  
  metrics_pland    <- landscapemetrics::lsm_c_pland(clip2)
  metrics_pland$id <- buffer$ID
  
  return(metrics_pland)
}

# Calculate proportion of missing cells

calc_prop_miss <- function(l, buffer) {
  
  clip1 <- raster::crop(l, extent(buffer))
  clip2 <- raster::rasterize(buffer, clip1, mask = TRUE)

  if(raster::cellStats(clip2 > 0, stat = "sum") == 0) {
    
    # if there are no populated cells
    metric_calc <- NA
    prop_missing <- data.frame(
      "id" = buffer$ID,
      "prop_missing" = 1)
    clip2 <- NA
    
  } else{
    
    # create mask
    clip3 <- rasterize(buffer, setValues(clip1, 1), mask = TRUE)
    
    # calculate metrics 
    prop_missing <- 1 - cellStats(!is.na(clip2), sum) / cellStats(clip3, sum)
    prop_missing <- data.frame(
      "id" = buffer$ID,
      "prop_missing" = prop_missing)
  }
  
 return(prop_missing)
  
}


calc_pland <- function(l, buffer) {
  
  clip1 <- raster::crop(l, extent(buffer))
  clip2 <- raster::rasterize(buffer, clip1, mask = TRUE)
  
  f <- freq(clip2)
  f <- data.frame(f)
  f <- subset(f, value %in% c(1:7))
  f$area <- f$count * 0.6 

  p <- data.frame(f, p = f[,2]/sum(f[, 2]))
  
  return(metrics_pland)
}
  
  