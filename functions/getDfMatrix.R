#updated on 26.08.2017
getDfMatrix <- function(dfMatrix, fdr_threshold = 0.1, lfc_threshold = log2(2)){
	colnames(dfMatrix) <- c("logFC","padj")
	dfMatrix[is.na(dfMatrix$logFC),"logFC"] <- 0.0 #converts all logFC NA values to 0.0, no expression!
	dfMatrix[is.na(dfMatrix$padj),"padj"] <- 1.0   #converts all padj NA values to 1.0, not significant!
	dfMatrix[dfMatrix$padj >= 0.1, "padj"] <- 1.0
	dfMatrix[dfMatrix$padj < 0.1, "padj"] <- 0.0
	dfMatrix[dfMatrix$padj >= fdr_threshold,"logFC"] <- 0.0
	dfMatrix[dfMatrix$logFC <= -lfc_threshold & dfMatrix$padj < fdr_threshold,"logFC"] <- -1.0
	dfMatrix[dfMatrix$logFC >= lfc_threshold & dfMatrix$padj < fdr_threshold,"logFC"] <- 1.0
	dfMatrix[abs(dfMatrix$logFC) < lfc_threshold,"logFC"] <- 0.0
	dfMatrix$padj <- NULL
	return(dfMatrix)
}