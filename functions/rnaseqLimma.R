#updated on 20.10.2017
rnaseqLimma <- function(counts = NULL, design = NULL, contrast.matrix = NULL){	
	# the rnaseqLimma functions requires 4 inputs:
	# (1) counts  : datafame with only raw reads and rownames
	# (2) design  : the design matrix 
	# (3) MatrixNames : The name of the contrast matrix for the output, A-B as A_vs_B!
	# (4) contrast.matrix  : The contrat matrix to determine the differentially expressed genes!
	# for details please conduct the limma user guide!
	dgeList <- DGEList(count = counts)
	dgeList <- calcNormFactors(dgeList)
	normVoom <- voom(dgeList, design, plot = TRUE, normalize = "quantile")
	fit <- lmFit(normVoom,design)
	fit <- contrasts.fit(fit, contrast.matrix)
	fit <- eBayes(fit)
	return(fit)
}