#updated on 26.08.2017
plotLogs <- function(orgDfr = NULL, fdr_threshold = 0.1, output_folder = NULL, legendpos = "bottomleft"){
	print("Saving The Logs vs DEG!...")
	subDfr <- orgDfr[,seq(1,ncol(orgDfr),2)]
	for(i in 1:ncol(subDfr)){subDfr[,i] <- 0}
	logs = seq(1.4,2.0,0.01)
	outDfr <- as.data.frame(matrix(nrow = length(logs), ncol = ncol(subDfr)*3))
	colnames(outDfr) <- paste(rep(colnames(subDfr),  each = 3),c("(-1)","(0)","(+1)"),sep = "_")
	rownames(outDfr) = as.character(logs)
	for(i in 1:length(logs)){
		lfc = logs[i]
		indexes = seq(2,ncol(orgDfr),2)
		hh = 1
		for(j in indexes){
			tmpTb = c(0,0,0)
			class(tmpTb) <- "table"
			names(tmpTb) <- c("-1","0","1")
			subDfr[,which(indexes == j)] = getDfMatrix(orgDfr[,(j-1):j],lfc_threshold = log2(lfc),fdr_threshold = 0.1)
			tb = table(subDfr[,which(indexes == j)])
			tmpTb[c(which(names(tmpTb) %in% names(tb)))] = tb
			tb = tmpTb
			outDfr[i,hh] = tb[[1]]
			hh = hh + 1
			outDfr[i,hh] = tb[[2]]
			hh = hh + 1
			outDfr[i,hh] = tb[[3]]
			hh = hh + 1
		}
	}
	outDfr <<- outDfr[,-(seq(2,ncol(outDfr),3))]
	#write.table(outDfr, file = "logs.tab",sep = "\t")
	cols = c("green","red")
	pchs = rep(seq(1,ncol(outDfr),2), each = 2)
	rownames(outDfr) = as.character(logs)
	legends = unique(gsub("\\..*","",colnames(outDfr)))
	legendCols = rep("black", length(legends))
	#legends = c(legends,"DOWN","UP")
	#legendCols = c(legendCols,cols)
	pdf(output_folder)
	max = max(outDfr)
	print(max)
	print(rowSums(abs(outDfr)))
	plot(logs,outDfr[,1], type = "l", pch = 1, cex = 1, xlab = "log2Fold", ylab = "Number of DE Genes", col = cols[1], xlim = c(0,2.1), ylim = c(0, max))
	for(i in 2:ncol(outDfr)){
		lines(logs, outDfr[,i], pch = pchs[i], cex = 1, col = cols[(2-(i%%2))])
	}
	legend(legendpos,bty = "n", text.width = 1, xjust=1, yjust=1, legend= legends, pch = c(unique(pchs),20,20), col= legendCols, cex = 0.8)
	legend("topright",bty = "n", text.width = 1, xjust=1, yjust=1, legend = c("DOWN","UP"), pch = 20, col= cols, cex = 0.8)
	dev.off()
}