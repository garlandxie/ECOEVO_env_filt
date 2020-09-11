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
  metrics_pland$id <- buffer$Site_ID
  
  return(metrics_pland)
}