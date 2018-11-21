######################################################################
# Extracting breeding range envrionmental data from WorldClim data   #
######################################################################

    
             
presab <- function(path, grid_size=1, crs=NULL, pres=1, orig=1, season=c(1,2), clip=TRUE){
	# Libraries
	require(sf)
	require(fasterize)
	require(dplyr)
	require(raster)
	require(maptools)
	data("wrld_simpl", package = 'maptools')
	world <- st_as_sf(wrld_simpl)
	
	
	# Create an empty raster. The default crs string creates a 1 degree equal area grid based on the WGS84 ellipsoid model. Grid cells are approximately 100km2. Other crs strings can be defined by the user.
	if (is.null(crs)) {crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"} 
	null_rast <- raster(crs=crs)
	extent(null_rast) <- c(-180, 180, -60, 90)
	res(null_rast) <- 1
	
	# Optionally rescale the grid. 
	if (grid_size > 1) {null_rast <- aggregate(null_rast, fact=grid_size)}
	if (grid_size < 1) {null_rast <- disaggregate(null_rast, fact=1/grid_size)}
	
	

	# Get list of shapefiles. It is assumed that all shapefiles (and only shapefiles) are in the path
	files <- list.files(paste(path), pattern=".shp")
	
	# Create empty presence absence matrix
	pres_ab <- matrix(NA, ncol=length(files), nrow=null_rast@ncols*null_rast@nrows)
	files <- gsub(".shp", "", files, fixed=T)
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
   
   worldrast<-fasterize(world, null_rast, field=NULL)
   crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
   crs(worldrast) <- crs    

   pres_ab[is.na(worldrast@data@values)] <- NA
   
   }
  return(list(pres_ab, crs, grid_size))
}





############# Extend to allow summaries of multiple traits e.g. sum of variance
# Species attribute to grid
sp_trait <- function(presab_dat, trait, cell_function="median", plot_map=TRUE, clip=TRUE, clip_by_richness=0, scale="eq_number", nbins=50){
	# Libraries
	require(sf)
	require(fasterize)
	require(dplyr)
	require(raster)
	require(dispRity)
	require(letsR)
	require(ggplot2)

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
			null_rast2@data@min <- min(null_rast2@data@values, na.rm=TRUE)
			null_rast2@data@max <- max(null_rast2@data@values, na.rm=TRUE)
			
			if(clip){
				data(temp)
				crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
				crs(temp) <- crs    
				temp <- aggregate(temp, fact=6)
				if (grid_size > 1) {temp <- aggregate(temp, fact=grid_size)}
				if (grid_size < 1) {temp <- disaggregate(temp, fact=1/grid_size)}
				null_rast2@data@values[is.na(temp@data@values)] <- NA
				
				if(clip_by_richness>0){
					rich <- sp_rich(birds_presab, plot_map=FALSE)
					null_rast2@data@values[which(rich[[1]]@data@values<clip_by_richness)] <- NA
				}
   data.values <- null_rast2@data@values
   }
   
	if(scale=="eq_number") { cuts <- cut_number(null_rast2@data@values, n = nbins)
		null_rast2@data@values <- as.numeric(cuts)}
	if(scale=="eq_width") { cuts <- cut_interval(null_rast2@data@values, n = nbins)
		null_rast2@data@values <- as.numeric(cuts)}
	if(scale=="raw") { null_rast2@data@values <- null_rast2@data@values; cuts <- NA}
			if (plot_map==TRUE) {plot(null_rast2)}			
				if (plot_map==TRUE) {plot(null_rast2)}
				return(list(null_rast2, cuts, data.values))
				}

	
	
# Species richness
sp_rich <- function(presab_dat, plot_map=TRUE, clip=TRUE, scale="eq_number", nbins=50){
	# Libraries
	require(sf)
	require(fasterize)
	require(dplyr)
	require(raster)
	require(letsR)
	
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
   
   data(temp)
   crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
   crs(temp) <- crs    
   temp <- aggregate(temp, fact=6)
   	if (grid_size > 1) {temp <- aggregate(temp, fact=grid_size)}
	if (grid_size < 1) {temp <- disaggregate(temp, fact=1/grid_size)}
   null_rast2@data@values[is.na(temp@data@values)] <- NA
   null_rast2@data@values[null_rast2@data@values==0] <- NA
   data.values <- null_rast2@data@values
   }
   
	if(scale=="eq_number") { cuts <- cut_number(null_rast2@data@values, n = nbins)
		null_rast2@data@values <- as.numeric(cuts)}
	if(scale=="eq_width") { cuts <- cut_interval(null_rast2@data@values, n = nbins)
		null_rast2@data@values <- as.numeric(cuts)}
	if(scale=="raw") { null_rast2@data@values <- null_rast2@data@values; cuts <- NA}
			if (plot_map==TRUE) {plot(null_rast2)}
	return(list(null_rast2, cuts, data.values))
	}




color.legend2 <- function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
    gradient = "x", border, ...) 
{
    oldcex <- par("cex")
    par(xpd = TRUE, cex = cex)
    gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), 
        gradient = gradient, border=border)
    if (gradient == "x") {
        xsqueeze <- (xr - xl)/(2 * length(rect.col))
        textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
        if (match(align, "rb", 0)) {
            texty <- yb - 0.2 * strheight("O")
            textadj <- c(0.5, 1)
        }
        else {
            texty <- yt + 0.2 * strheight("O")
            textadj <- c(0.5, 0)
        }
    }
    else {
        ysqueeze <- (yt - yb)/(2 * length(rect.col))
        texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
        if (match(align, "rb", 0)) {
            textx <- xr + 0.2 * strwidth("O")
            textadj <- c(0, 0.5)
        }
        else {
            textx <- xl - 0.2 * strwidth("O")
            textadj <- c(1, 0.5)
        }
    }
    text(textx, texty, labels = legend, adj = textadj, ...)
    par(xpd = FALSE, cex = oldcex)
}

