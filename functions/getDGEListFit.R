#updated 26.08.2017
getDGEListFit <- function(dgeLists,sampleNameOrdered, conditions) {
	#rownames(dgeLists) <- dgeLists$GeneID
	#dgeLists$GeneID <- NULL
	dgeList <- DGEList(counts = dgeLists, genes = rownames(dgeLists), group = conditions)
	#dgeList$samples$lib.size <- colSums(dgeList$counts)
	dgeList <- calcNormFactors(dgeList)
	design <- model.matrix(~0+conditions, data = dgeList$samples)
	rownames(design) <- colnames(dgeList)
	colnames(design) <- sampleNameOrdered
	dgeList <- estimateDisp(dgeList, design)
	fit <- glmFit(dgeList, design)
	return(fit)
}