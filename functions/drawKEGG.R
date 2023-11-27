drawKEGG <- function(dfMatrix = NULL, keggids = NULL, output_folder = NULL, gene.idtype = NULL, keggNative = TRUE){
	setwd(output_folder)
	path_back <- paste(rep("..",(str_count(output_folder,"/")+1)),collapse = "/")
	org3L <- orgInfo$org3L
	keggDir <- paste0(mainDir,"/keggPathDir/",paste0(org3L,"PNGkgml"),"/")
	if(!("keggPathDir" %in% dir(mainDir))) dir.create(paste0(mainDir,"/keggPathDir/"),showWarnings = FALSE)
	if(!(paste0(org3L,"PNGkgml") %in% dir(paste0(mainDir,"/keggPathDir/")))){
		getPathPNGkgml(orgId = orgInfo$orgId, org3L = orgInfo$org3L, outputFolder = keggDir)
	}
	keggids <- keggids[!(keggids$KEGGID %in% "01100"),]
	keggids <- keggids[keggids$Count > 0,]
	keggids <- keggids[keggids$Pvalue < 0.05,]
	keggids$Term <- gsub("[/:; ]","",keggids$Term)
	for(i in 1:nrow(keggids)){
	  # use try, as for some pathway
	  try(pv.out <- pathview(gene.data = dfMatrixDfr, pathway.id = keggids$KEGGID[i], kegg.dir = keggDir, species = org3L, gene.idtype = gene.idtype, out.suffix = gsub(" ","_",keggids$Term[i]), kegg.native = keggNative))
	}
	setwd(path_back)
}
