######################################################################
# Create presence absence matrix using fasterize. Note fasterize counts presence 
# only if polygon overlaps the cell centroid so can miss some presences.
# The presabBL function samples at the 0.5 degree scale with the option to aggregate at courser
# scales. To retain speed it does this on the full presence-absence matrix at the end the run.
# Finer scales are not possible with this approach because of the memory cost of storing a massive matrix.
# 
# TO DO
# Finer scales should be possible with finer grids applied to each shapefile followed by aggregation to the PA matrix.
# Keep an eye on updates to fasterize for potential addition of getCover (would simplify these resceling issues)
######################################################################

presabBL <- function(path, species=NULL, grid_size=1, crs=NULL, pres=1, orig=1, season=c(1,2), clip=FALSE, longlat_extent=c(-180, 180, -90, 90)){
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
  
  null_rast <- disaggregate(null_rast, fact=1/0.5)
  
  # Get list of shapefiles. It is assumed that all shapefiles (and only shapefiles) are in the path
  files <- list.files(paste(path), pattern=".shp")
  files <- gsub(".shp", "", files, fixed=T)
  
  if (!is.null(species)){files <- intersect(species, files)}
  
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
  
  if(grid_size %in% c(0.5,1,2,4) == FALSE) {grid_size <- 1; warning("Grid size must be 0.5, 1, 2, or 4. Setting grid_size to 1.")}
  
  if (grid_size==0.5) { return(list(pres_ab, crs, grid_size))  }
  
  if (grid_size %in% c(1,2,4)) {
    
    x<-as.numeric(rep(gl(360/grid_size, grid_size*2), grid_size*2))
    y <- c()
    for (i in 1:(180/grid_size)){y <- c(y, x+(360)*(i-1))}
    pres_ab <- rowsum(pres_ab, as.integer(y), na.rm=TRUE)
    pres_ab[pres_ab>0] <- 1
    pres_ab[pres_ab==0] <- NA
    return(list(pres_ab, crs, grid_size))
  }
}