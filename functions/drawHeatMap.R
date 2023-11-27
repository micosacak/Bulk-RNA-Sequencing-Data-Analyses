#updated on 20.10.2017
drawHeatMap <- function(theMatrix = NULL, output_folder = NULL, main = "heatMap"){
	pdf(paste0(output_folder,"/heatMap.pdf"))
  showLabRow = TRUE
  if(nrow(theMatrix) >= 60) showLabRow = NA
	nCol = ncol(theMatrix)
	mycols <- brewer.pal(8, "Dark2")[1:length(unique(nCol))]			
	heatmap.3(theMatrix, 
		key = T, keysize=1, density.info = c("none"), denscol="yellow", trace="none",
		col=colorpanel(100, "green","yellow", "red"), margin = c(length(nCol)*10, 
		length(nCol)*25), main= main, cexRow = 0.6, cexCol = 0.5, labRow = showLabRow)	
	dev.off()
}