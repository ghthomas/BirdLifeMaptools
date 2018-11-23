plot_bird_raster <- function(raster_data=NULL, scale="raw", nbins=20, col=viridis(20), cex=1.5, legend_pos=c(-80,-80,80, -75), world=TRUE){
  require(plotrix)
  require(ggplot2)
  require(viridis)
  
  if(scale=="eq_number") { cuts <- cut_number(raster_data@data@values, n = nbins)
  raster_data@data@values <- as.numeric(cuts)}
  if(scale=="eq_width") { cuts <- cut_interval(raster_data@data@values, n = nbins)
  raster_data@data@values <- as.numeric(cuts)}
  if(scale=="raw") {  cuts <- cut_interval(raster_data@data@values, n = length(na.omit(unique(raster_data@data@values))))
  raster_data@data@values <- as.numeric(cuts)}
  
  cuts_vec <- sort(na.omit(unique(cuts)))
  idx <- c(round(quantile(1:length(cuts_vec), c(0,0.2,0.4,0.6,0.8,1))))
  brks <- as.character(cuts_vec[idx])
  brks <- gsub("\\[|\\]|\\)|\\(", "", brks)
  brks <- as.numeric(unlist(strsplit(brks, ",")))[c(1,3,5,7,9,12)]
 # brks[1] <- raster_data@data@min
#  brks[6] <- raster_data@data@max
  
  if (world==TRUE) {
    require(maptools)
    data("wrld_simpl", package = 'maptools')
    plot(wrld_simpl, col="light grey", border=NA, bg="white")
    par(bty="n")
    plot(raster_data, add=T, legend=F, col=col, axes=FALSE)
    }
  if (world==FALSE) {par(bty="n")
    plot(raster_data[[1]], add=F, legend=F, col=col, axes=FALSE)
  }

  
  color.legend2(legend_pos[1], legend_pos[2], legend_pos[3], legend_pos[4], brks, rect.col=col, bty="n", align="rb", cex=1.5, border=NA)
  
}
  