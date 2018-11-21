sp_rich <- function(presab_dat, plot_map=FALSE, clip=TRUE){
  # Libraries
  require(sf)
  require(fasterize)
  require(dplyr)
  require(raster)

  # Create an empty raster
  null_rast <- raster(crs=crs)
  extent(null_rast) <- c(-180, 180, -60, 90)
  res(null_rast) <- 1
  
  grid_size <- presab_dat[[3]]
  
  # Optionally rescale the grid. 
  if (grid_size > 1) {null_rast <- aggregate(null_rast, fact=grid_size)}
  if (grid_size < 1) {null_rast <- disaggregate(null_rast, fact=1/grid_size)}
  
  null_rast2 <- null_rast
  null_rast2@data@inmemory <- TRUE
  null_rast2@data@values <- rowSums(presab_dat[[1]], na.rm=TRUE)
  null_rast2@data@min <- min(null_rast2@data@values)
  null_rast2@data@max <- max(null_rast2@data@values)

  if(clip){
    data("wrld_simpl", package = 'maptools')
    world <- st_as_sf(wrld_simpl)
    worldrast<-fasterize(world, null_rast, field=NULL)
    crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(worldrast) <- crs
    null_rast2@data@values[is.na(worldrast@data@values)] <- NA
  }
  
  return(null_rast2)
}