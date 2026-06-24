# extract_mean_raster_val.R

# extract mean raster value for each admin

library(terra)          # for rasters and vector data
library(dplyr)          # for tidy data manipulation
library(exactextractr)  # very efficient zonal stats (handles partial coverage)



# 3. Function to extract mean pixel value for one raster
extract_mean = function(raster_path, admin_shape, admin_colname, layer=NA) {
  raster = rast(raster_path)
  # if the raster has multiple layers, use the specified one
  if(grepl('grd', raster_path)){
    # Select the appropriate raster layer
    raster = raster[[layer]]
  } else if (!grepl('tif', raster_path)) warning('Did not recognize raster file type. Expected .tif or .grd')

  # Make sure CRS matches shapefile
  if (!(crs(raster) == crs(admin_shape))) {
    admin_shape = project(admin_shape, crs(raster))
  }
  
  if (!identical(crs(raster), crs(admin_shape))) {
    admin_shape = project(admin_shape, crs(raster))
  }
  
  # Extract mean values by county
  vals = terra::extract(raster, admin_shape, 'mean', na.rm=TRUE)
  
  # Combine with county ID and year
  out = data.frame(
    county_id = admin_shape[[admin_colname]],
    mean_val  = vals,
    year      = gsub("\\D", "", basename(raster_path))
  )
  return(out)
}
