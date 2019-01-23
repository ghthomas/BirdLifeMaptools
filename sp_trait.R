sp_trait <- function(presab_dat, trait, cell_function="median", range_weights=NULL, plot_map=TRUE, clip=TRUE, clip_by_richness=0){
  # Libraries
  require(sf)
  require(fasterize)
  require(dplyr)
  require(raster)
  require(dispRity)
  require(ggplot2)
  
  if (is.null(rownames(trait))) {stop("Trait matrix must be supplied with rownames.")}
  
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
  
  
  if(!is.matrix(trait)) {stop("Trait data must be supplied as a matrix")}
  
  if(cell_function =="sumvar" & dim(trait)[2]<2) {stop("At least two traits are required to calculate the sum of variance")}
  if(cell_function =="centroid_dist" & dim(trait)[2]<2) {stop("At least two traits are required to calculate the mean centroid distance")}
  if(cell_function =="convhull" & dim(trait)[2]<2) {stop("At least two traits are required to calculate the convex hull volume")}
  
  tr_val <- rep(NA, dim(presab_dat[[1]])[1])
  for (i in 1:dim(presab_dat[[1]])[1]) {
    tmp_val <- NA
    sp_index <- !is.na(presab_dat[[1]][i,])
    tmp_nms <- colnames(presab_dat[[1]])[sp_index]
    idx <- which(rownames(trait) %in% tmp_nms) 
    
    if (length(idx)>0){
      trait_vec <- matrix(trait[idx,], ncol=dim(trait)[2])
      
      if (cell_function =="median") {tmp_val <- median(trait_vec, na.rm=TRUE)}
      if (cell_function =="mean") {tmp_val <- mean(trait_vec, na.rm=TRUE)}
      if (cell_function =="variance") {tmp_val <- var(trait_vec, na.rm=TRUE)}
      
      if(dim(trait_vec)[1]>1){
        if (cell_function =="sumvar") {tmp_val <- sum(apply(trait_vec, 2, var, na.rm=TRUE))}
        if (cell_function =="centroid_dist") {tmp_val <- mean(dispRity(trait_vec, metric=centroids)$disparity[[1]][[1]])}
        if(dim(trait_vec)[2]<dim(trait_vec)[1]){
          if (cell_function =="convhull") {tmp_val <- dispRity(trait_vec, metric= convhull.volume)$disparity[[1]][[1]]}
        }
      }
      
      tr_val[i] <- tmp_val
    } else {tr_val[i] <- NA}
  }
  
  null_rast2@data@values <- tr_val
  null_rast2@data@min <- min(tr_val, na.rm=TRUE)
  null_rast2@data@max <- max(tr_val, na.rm=TRUE)
  
  if(clip){
    data("wrld_simpl", package = 'maptools')
    world <- st_as_sf(wrld_simpl)
    worldrast<-fasterize(world, null_rast, field=NULL)
    crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    crs(worldrast) <- crs
    null_rast2@data@values[is.na(worldrast@data@values)] <- NA
  }
  
    
    if(clip_by_richness>0){
      rich <- sp_rich(presab_dat, plot_map=FALSE)
      null_rast2@data@values[which(rich@data@values<clip_by_richness)] <- NA
    }
  

  if (plot_map==TRUE) {plot(null_rast2)}			
  return(null_rast2)
}
