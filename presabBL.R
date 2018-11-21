presabBL <- function(path, species=NULL, grid_size=1, crs=NULL, pres=1, orig=1, season=c(1,2), clip=TRUE, longlat_extent=c(-180, 180, -60, 90)){
  # Libraries
  require(sf)
  require(fasterize)
  require(dplyr)
  require(raster)
  require(maptools)
  
  if (is.null(species)){ print("No species list provided, analysing all files in path.")}
  if (!is.null(species)){ print("Analysing species list.")}
  
    # Create an empty raster. The default crs string creates a 1 degree equal area grid based on the WGS84 ellipsoid model. Grid cells are approximately 100km2. Other crs strings can be defined by the user.
  if (is.null(crs)) {crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"} 
  null_rast <- raster(crs=crs)
  extent(null_rast) <- longlat_extent
  res(null_rast) <- 1
  
  # Optionally rescale the grid. 
  if (grid_size > 1) {null_rast <- aggregate(null_rast, fact=grid_size)}
  if (grid_size < 1) {null_rast <- disaggregate(null_rast, fact=1/grid_size)}
  
  # Get list of shapefiles. It is assumed that all shapefiles (and only shapefiles) are in the path
  files <- list.files(paste(path), pattern=".shp")
  files <- gsub(".shp", "", files, fixed=T)
  files <- intersect(species, files)
  
  if (length(files)!=length(species)) {"Some names in species list do not match names in the file path. Analysis will run but you may need to check your species list."}
  
  # Create empty presence absence matrix
  pres_ab <- matrix(NA, ncol=length(files), nrow=null_rast@ncols*null_rast@nrows)
  colnames(pres_ab) <- files
  
  # Now the important stuff
  for (i in 1:length(files)){
    fspp <- files[i]
    
    dat <- st_read(dsn=paste(path), layer=fspp, quiet=TRUE)
    dat <- dat %>% filter(presence %in% pres, origin %in% orig, seasonal %in% season)
    if (nrow(dat) > 0) {
      dat <- st_transform(dat, crs=crs)
      sp_rast <- fasterize(dat, null_rast, field=NULL)
      idx <- !is.na(sp_rast@data@values)
      pres_ab[idx, i] <- rep(1, sum(idx))
    } else {newval <- NA}
    
    cat(i, "of", length(files), "\r") 
    flush.console()
  }
  
  if(clip){
    data("wrld_simpl", package = 'maptools')
    world <- st_as_sf(wrld_simpl)
    worldrast<-fasterize(world, null_rast, field=NULL)
    crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(worldrast) <- crs
    pres_ab[is.na(worldrast@data@values)] <- NA
  }
  return(list(pres_ab, crs, grid_size))
}