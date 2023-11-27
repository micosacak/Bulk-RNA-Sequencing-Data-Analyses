#updated on 20.10.2017
plotSampleDists <- function(rld, output_folder, sub_conditions){
	sampleNumber = length(unique(sub_conditions))
	if(sampleNumber >= 3) {
		mycols <- brewer.pal(sampleNumber, "Set2")
	}else{
		mycols = c("darkgreen","darkorange")
	}
	
	sampleDists <- as.matrix(dist(t(assay(rld))))
	pdf(paste(output_folder,"heatmap_sample2sample_distances.pdf", sep = "/"))
	heatmap.3(sampleDists, 
	key = T, keysize=1, density.info = c("none"), trace="none",
	col = colorpanel(200, "black", "white"), 
	RowSideColors = t(as.matrix(mycols[as.integer(as.factor(sub_conditions))])), 
	ColSideColors = as.matrix(mycols[as.integer(as.factor(sub_conditions))]),
	margins = c(13,14))
	#legend("topright",legend = unique(sub_conditions),col = mycols[as.integer(as.factor(unique(sub_conditions)))], pch = 20, cex =0.5)
	plot.new()
	legend("center",legend = unique(sub_conditions),col = mycols[as.integer(as.factor(unique(sub_conditions)))], pch = 20, cex =1)
	
	dev.off()
}
