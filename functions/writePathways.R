#updated on 20.10.2017
writePathways <- function(pathIds = NULL, orgDb = NULL, output_folder = NULL, resultDfr = NULL, geneIDtype = "ENSEMBL", prefix = ""){
	keyid2PATH <- select(orgDb, keys = keys(orgDb,keytype = geneIDtype), column = c("SYMBOL","GENENAME","PATH"), keytype = geneIDtype)
	colnames(keyid2PATH) <- c("geneIds","SYMBOL","GENENAME","PATH")
	pathNames <- toTable(KEGGPATHID2NAME)
	keyid2PATH <- merge(keyid2PATH,pathNames, by.x = "PATH", by.y = "path_id", all.x = TRUE, all.y = FALSE)
	keyid2PATH <- keyid2PATH[!is.na(keyid2PATH$PATH),]
	printouts <- keyid2PATH[keyid2PATH$PATH %in% pathIds,]
	printouts <- printouts[!duplicated(printouts$PATH),]
	lt = vector("list",length = length(pathIds))
	names(lt) <- paste(orgInfo$org3L,pathIds, sep = "")
	for(i in 1:length(names(lt))){
		lt[[i]] <- data.frame()
	}
	for(pathId in pathIds){
		geneIds <- keyid2PATH$geneIds[keyid2PATH$PATH == pathId]
		outFile <- paste0(orgInfo$org3L,pathId)
		tmpRes <- resultDfr[rownames(resultDfr) %in% geneIds,]
		tmpRes <- merge(keyid2PATH, tmpRes, all.x = FALSE, all.y = TRUE, by.x = "geneIds", by.y = "row.names")
		tmpRes <- tmpRes[tmpRes$PATH == pathId,]
		head(tmpRes)
		if(nrow(tmpRes) > 0){
			lt[outFile][[1]] <- tmpRes
		}else{
			names(lt[which(names(lt) == outFile)]) <- strtrim(names(lt[which(names(lt) == outFile)]),3)
		}
	}
	chooseNames = c()
	for(nnn in names(lt)) if(nrow(lt[nnn][[1]]) > 0) chooseNames = c(chooseNames,nnn)
  lt = lt[names(lt) %in% chooseNames]
	write.xlsx(lt, file = paste0(output_folder,"/",prefix,"_keggGenes.xlsx"))
}
