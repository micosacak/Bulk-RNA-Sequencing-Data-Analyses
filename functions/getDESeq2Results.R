# 27.10.2017
getDESeq2Results <- function(aCo = NULL, dds = NULL, lfc_threshold = log2(2), fdr_threshold = 0.1){
  rs <- vector("list", length = 1)
	names(rs) <- c("deseq")
	for(i in 1:3){
		rs[[i]] <- vector("list", length = length(aCo))
		names(rs[[i]]) <- aCo
	}
	for(i in 1:length(aCo)){
		output_folder <- aCo[i]
		conditionSplit <- strsplit(output_folder,"_vs_")
		trt_sp = conditionSplit[[1]][1]
		unt_sp = conditionSplit[[1]][2]
		rs$deseq[[i]] <- as.data.frame(results(dds, contrast = c("condition",trt_sp, unt_sp)))
	}
	dfr <- vector("list", length = 1)
	names(dfr) <- c("deseq")
	dfr[[1]] <- vector("list", length = 2)
	names(dfr[[1]]) <- c("outDfr","pathDfr")
	tmpDfr <- vector("list", length = length(aCo))
	names(tmpDfr) <- aCo
	for(j in 1:length(aCo)){
		tmpDfr[[j]] <- rs$deseq[[j]]
		colnames(tmpDfr[[j]]) <- paste(aCo[j], colnames(tmpDfr[[j]]), sep = ".")
	}
	#deseq2  data!
	dfr$deseq$outDfr <- mergeAllDfr(tmpDfr)
	dfr$deseq$pathDfr <- as.data.frame(matrix(0,nrow = nrow(dfr$deseq$outDfr), ncol = length(aCo)))
	rownames(dfr$deseq$pathDfr) <- rownames(dfr$deseq$outDfr)
	colnames(dfr$deseq$pathDfr) <- aCo
	cc = 1
	for(i in seq(2,6*length(aCo),6)){
		dfr$deseq$pathDfr[,cc] <- getDfMatrix(dfr$deseq$outDfr[,c(i,i+4)], fdr_threshold = 0.1, lfc_threshold = lfc_threshold)
		cc <- cc + 1
	}
	return(list("rs" = rs, "dfr" = dfr))
}
