#updated on 20.10.2017
runDESeq2dds <- function(trt_sp,unt_sp, ddsFull){
	dds <- ddsFull[,colData(ddsFull)$condition %in% c(trt_sp,unt_sp)]
	dds$condition <- droplevels(dds$condition)
	dds$condition <- relevel(dds$condition,unt_sp)
	dds <- DESeq(dds, parallel = T)
	return(dds)
}