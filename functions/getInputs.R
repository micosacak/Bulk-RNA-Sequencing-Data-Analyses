getInputs <- function(rawCounts = NULL, comparisons = NULL){
	all_conditions = gsub("\\_R[0-9]","",colnames(rawCounts))
	conditions = factor(as.integer(factor(all_conditions)))
	design <- model.matrix(~0 + conditions)
	colnames(design) <- levels(factor(all_conditions, levels = unique(all_conditions)))
	#print(colnames(design))
	comps = c()
	for(i in seq(2,length(comparisons),2)){
		comps = c(paste0(colnames(design)[comparisons[i-1]],"-",colnames(design)[comparisons[i]]),comps)
	}
	contrast.matrix <- makeContrasts(contrasts = comps, levels= design)
	return(list("rawCounts"=rawCounts, "all_conditions" = all_conditions, "conditions" = conditions, "design" = design, "contrast.matrix" = contrast.matrix, "sampleNameOrdered" = colnames(design)))
}