#updated on 20.10.2017
writeGOs <- function(goDfr = NULL, orgDb = NULL, output_folder = NULL, resultDfr = NULL, geneIDtype = "ENSEMBL", goType = NULL){
	orgDb <- eval(parse(text = orgDb))
	keyid2GO <- select(orgDb, keys = keys(orgDb,keytype = geneIDtype), column = c("GO"), keytype = geneIDtype)
	keyid2GO <- keyid2GO[,-c(3:4)]
	head(keyid2GO)
	colnames(keyid2GO) <- c("geneIds","GO")
	tmpVector = as.character(unname(unlist(goDfr[1,])))
	colnames(goDfr)[which(startsWith(tmpVector,"GO:"))] <- "GO"
	colnames(goDfr)[which(tolower(colnames(goDfr)) == "term")] <- "Term"
	if(nrow(goDfr) != 0){
		lt <- vector("list",length = length(goDfr$GO))
		names(lt) <- gsub("[/:]","",goDfr$GO)
		for(i in 1:length(names(lt))){
			lt[[i]] <- data.frame()
		}
		for(j in 1:nrow(goDfr)){
			goId <- goDfr$GO[j]
			geneIds <- keyid2GO$geneIds[keyid2GO$GO == goId]		
			tmpRes <- resultDfr[rownames(resultDfr) %in% geneIds,]
			tmpGoId <- gsub(":","",goId)
  		if(nrow(tmpRes) > 0){
  			prefix <- goDfr$Term[goDfr$GO == goId]
  			outFile <- paste0(tmpGoId,"_",goDfr$Term[j])
  			outFile  <- gsub("[/ :]","", outFile) #replace any special character with ""
  			outFile <- strtrim(outFile,30)        #trimming might be better to prevent of writing into files!	
  			tmpRes <- merge(keyid2GO, tmpRes, all.x = FALSE, all.y = TRUE, by.x = "geneIds", by.y = "row.names")
  			tmpRes <- tmpRes[tmpRes$GO == goId,]
  			tmpRes[,"Term"] <- rep(goDfr$Term[j], each = nrow(tmpRes))
  			lt[tmpGoId][[1]] <- tmpRes
  			names(lt)[which(names(lt) == tmpGoId)] <- outFile
  		}else{
  			lt[which(names(lt) == tmpGoId)] <- NULL
  		}
	}
	lt <- lt[nchar(names(lt)) != 9]
	if(length(lt) >= 1) write.xlsx(lt, file = paste0(output_folder,"/",paste0(goType,"_genes.xlsx")))
	}
}
