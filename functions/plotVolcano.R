#updated on 20.10.2017
plotVolcano <- function(Df, topDeg = NA, plot3colors = FALSE, lfc_threshold = log2(1.5), fdr_threshold = 0.1, main = NULL, legendpos = "topleft", labelsig = FALSE, textcx = 0.5, xNames = "logFC", yNames = "-log10(p-value)", fdrNames = "FDR", ...) {
	main <- unlist(strsplit(main,split = "/"))
	if(length(main) > 1) main <- main[length(main)]
	colnames(Df) <- c("xVals","yVals","FDR")
	maxLim <- max(abs(Df$xVals[!is.na(Df$xVals)])) + 1
	upR <- nrow(subset(Df,(xVals) >= lfc_threshold & FDR < fdr_threshold))
	downR <- nrow(subset(Df,(xVals) <= -lfc_threshold & FDR < fdr_threshold))
	if(length(topDeg) <= 1){
		topDeg <- subset(Df,abs(xVals) >= lfc_threshold & FDR < fdr_threshold)
		topDeg <- topDeg[order(topDeg$FDR, decreasing = FALSE),]
		if(nrow(topDeg) < 20){
			print(c("nrow(topDeg) < ",20, " = ", nrow(topDeg)))
		}
		topDeg <- topDeg[1:20,]
	}else{
		topDeg <- subset(Df,rownames(Df) %in% topDeg)
	}
	if(plot3colors){
		with(Df, plot(xVals, -log10(yVals), pch = 20, cex = 0.5, xlim = c(-maxLim,maxLim), main = main, xlab = xNames, ylab = yNames, cex.main = 0.8, ...))
		with(subset(Df, FDR < fdr_threshold  & (xVals) >= lfc_threshold), points(xVals, -log10(yVals), pch = 20, cex = 0.5, col = "red", ...))
		with(subset(Df, FDR < fdr_threshold  & (xVals) <= -lfc_threshold), points(xVals, -log10(yVals), pch = 20, cex = 0.5, col = "green", ...))
		#with(subset(topDeg,abs(xVals) >= lfc_threshold & FDR < fdr_threshold), points(xVals, -log10(yVals), pch = 1, cex = 1, col = "black", ...))
		legend(legendpos,bty = "n", text.width = 1, xjust=1, yjust=1, legend=c("No Change", paste0("Up Regulated [",upR,"]"), paste0("Down Regulated [",downR,"]")), pch = 20, col=c("black","red","green"), cex = 0.8)
		#legend("top", bty = "n", text.width = 1, legend = c(paste0("up Regulated = ", upR), paste0("down Regulated = ",downR)), pch = 20, col=c("red","green"), cex = 0.5)
		print(c(upR,downR))
	}else{
		with(Df, plot(xVals, -log10(yVals), pch = 20, cex = 0.5, xlim = c(-maxLim,maxLim), main = main, xlab = xNames, ylab = yNames, ...))
		with(subset(Df, abs(xVals) >= lfc_threshold), points(xVals, -log10(yVals), pch = 20, cex = 0.5, col = "red", ...))  
		with(subset(Df, FDR < fdr_threshold), points(xVals, -log10(yVals), pch = 20, cex = 0.5, col = "pink", ...))
		with(subset(Df, FDR < fdr_threshold  & abs(xVals) >= lfc_threshold), points(xVals, -log10(yVals), pch = 20, cex = 0.5, col = "green", ...))
		legend(legendpos,bty = "n", text.width = 1, xjust=1, yjust=1, legend=c(paste(fdrNames,">", fdr_threshold, sep = ""), paste("|",xNames,"|>",round(lfc_threshold,2),sep=""), paste(fdrNames,"<",fdr_threshold ,sep=""), "Both"), pch = 20, col=c("black","red","pink","green"), cex = 0.5)
		#legend("bottomleft",bty = "n", text.width = 1, xjust=1, yjust=1, legend = c(paste0("up = ", upR), paste0("down = ",downR)), pch = 20, col=c("red","green"), cex = 0.5)
	}
}
