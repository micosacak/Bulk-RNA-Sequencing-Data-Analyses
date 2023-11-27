#updated on 20.10.2017
topGOanalysis <- function(uniGenes = NULL, output_folder = NULL, minNodeSize = 5, annotFunction = annFUN.org, ids = "ensembl", degDfr = NULL, ...){
	for(goTerm in c("BP","CC","MF")){	
		print(goTerm)
		sampleGOdata <- new(
		"topGOdata", 
		description = output_folder,
		ontology = goTerm,
		allGenes = uniGenes,
		geneSel = topDiffGenes,
		nodeSize = minNodeSize,
		annot = annotFunction,
		mapping = orgInfo$orgDb,
		ID = ids)
		print(c("Now analysing",goTerm," by topGO ..."))
		runTopGO(sampleGOdata, goTerm, output_folder, degDfr = degDfr)
	}
}