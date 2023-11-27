getGOframe <- function(geneIDtype, egGOdb, orgDb){
	if(geneIDtype == "ENTREZ"){
		dframe = toTable(egGOdb)
	}else if(geneIDtype == "ENSEMBL"){
		dframe <- AnnotationDbi::select(orgDb, keys = keys(orgDb, keytype = "ENSEMBL"), columns = c("GO"), keytype = "ENSEMBL")
		colnames(dframe) <- c("gene_id", "go_id", "Evidence", "Ontology")
		dframe <- dframe[!is.na(dframe$go_id),]
	}else{
		stop("gene id type must be ENSEMBL or ENTREZ")
	}
	return(dframe)
}