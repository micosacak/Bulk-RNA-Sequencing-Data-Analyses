#updated on 26.08.2017
plotMASmears <- function(Df, legendpos = "topleft", fdr_threshold = 0.1, lfc_threshold = log2(1.5), main = NULL, labelsig = TRUE, textcx = 1, xNames = "basemean",yNames = "logFC",fdrNames = "FDR",LoG = FALSE,... ){
	main <- unlist(strsplit(main,split = "/"))
	if(length(main) > 1) main <- main[length(main)]
	colnames(Df) <- c("xVals","yVals","FDR")
	if(LoG){
		Df$xVals <- log10(Df$xVals)
	}
	maxLim <- max(abs(Df$yVals[!is.na(Df$yVals)])) + 1
	with(Df, plot(xVals, yVals, pch = 20, cex= 0.5, ylim = c(-maxLim,maxLim), main = main, xlab = xNames, ylab = yNames, ...))
	with(subset(Df,(FDR < fdr_threshold & yVals >= lfc_threshold)), points(xVals,yVals, col = "red", pch = 20, cex = 0.5))
	with(subset(Df,(FDR < fdr_threshold & abs(yVals) < lfc_threshold)), points(xVals,yVals, col = "cyan", pch = 20, cex = 0.5))
	with(subset(Df,(FDR < fdr_threshold & yVals <= -lfc_threshold)), points(xVals,yVals, col = "green", pch = 20, cex = 0.5))
	#if (labelsig) {
	#	with(subset(Df, abs(yVals) >= lfc_threshold & FDR < sig_threshold), textxy(xVals, yVals, labs = rownames(Df), cex = textcx, col = 1000))
	#}
	legend(legendpos, cex = 0.5, bty = "n", xjust = 1, yjust = 1, legend = c(paste("FDR > ", fdr_threshold, sep=""), paste("FDR < ",fdr_threshold," & ", yNames," >= ",round(lfc_threshold,2), sep=""),paste("FDR < ",fdr_threshold, " & ", "|",yNames,"| < ", round(lfc_threshold,2), sep=""),  paste("FDR < ",fdr_threshold, " & ", yNames," <= -",round(lfc_threshold,2), sep="")), pch = 20, col=c("black","red","cyan","green"))
	abline( h = c(-lfc_threshold,0,lfc_threshold), col = "blue")
}