getPathPNGkgml <- function(orgId = "Dr",org3L = "dre",outputFolder = "."){
	#load required packages!
	library(paste0("org.",orgId,".eg.db"),character.only = TRUE)
	library("png")
	library("KEGGREST")
	library("XML")	
	pathids <- toTable(eval(parse(text = paste("org",orgId,"egPATH", sep = "."))))[,2]
	pathids <- unique(pathids)
	pathids <- paste0(org3L,pathids, collapse = ",")
	pathids <- unlist(strsplit(pathids,split = ","))
	dir.create(outputFolder, showWarnings = FALSE)
	noPaths <- c() # print keggids with no path kgml or image!
	for(keggid in pathids){
		try(pathPNG <- keggGet(keggid,"image"),silent = T)
		try(pathKGML <- keggGet(keggid,"kgml"), silent = T)
		if(length(pathPNG) != 0){
			keggXML <- xmlParse(pathKGML, asText = TRUE)
			saveXML(keggXML, file = paste0(outputFolder,"/",keggid,".xml"))
			writePNG(pathPNG, paste0(outputFolder,"/",keggid,".png"))
		}else{
			noPaths <- c(noPaths,keggid)
		}
	}
	print("Could not download image or kgml for the followings:\n")
	print(noPaths)
}
