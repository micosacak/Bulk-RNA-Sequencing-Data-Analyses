#updated on 20.10.2017
runTopGO <- function(sampleGOdata, goTerm, output_folder, degDfr = NULL){
	dir.create(output_folder ,showWarnings = FALSE)
	# algorithms and statistics that can be used together!
	#algoComp <- rbind(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
	#               c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
	#               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
	#               c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
	#               c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
	#               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
	# algorithms and statistics that have been tested!
	#algoComp <- 
	#	rbind(
	#	c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
	#	c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
	#	c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
	#	c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
	#	c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
	#	c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
	#	)
	# algorithms and statistics chosen to be used!	
	algoComp <- 
	rbind(
	c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
	c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
	c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
	c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
	c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
	c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
	)
	rownames(algoComp) <- c("classic", "elim", "weight", "weight01", "lea", "parentchild")
	colnames(algoComp) <- c("fisher", "z", "ks", "t", "globaltest", "category", "sum", "ks.ties")
	#algoComp
	if(sum(algoComp) == 0){
		stop("ERROR: no algorithms and test found to use!")
	}else{
		results <- vector("list", sum(algoComp))
	}
	h = 0
	for(i in 1:nrow(algoComp)){
		for(j in 1:ncol(algoComp)){
			if(algoComp[i,j]){
				h = h + 1
				names(results)[h] <- paste(rownames(algoComp)[i], colnames(algoComp)[j], sep = "")
				print(names(results)[h])
				results[h] <- runTest(sampleGOdata, algorithm = rownames(algoComp)[i], statistic = colnames(algoComp)[j], cutOff = 0.01)
			}
		}
	}
	allGO = usedGO(object = sampleGOdata) 
	# There are 48 possible comparisons with different statics and algorithms! User must define the allRes length using the following template!!!
	#allRes <- GenTable(sampleGOdata, 
	#				results[[1]], results[[2]], results[[3]], results[[4]], results[[5]],
	#				results[[6]], results[[7]], results[[8]], results[[9]], results[[10]],
	#				results[[11]], results[[12]], results[[13]], results[[14]], results[[15]],
	#				results[[16]], results[[17]], results[[18]],results[[19]], results[[20]], 
	#				results[[21]], results[[22]], results[[23]],results[[24]], results[[25]],
	#				results[[26]], results[[27]], results[[28]],results[[29]], results[[30]],
	#				results[[31]], results[[32]], results[[33]],results[[34]], results[[35]],
	#				results[[36]], results[[37]], results[[38]],results[[39]], results[[40]],
	#				results[[41]], results[[42]], results[[43]],results[[44]], results[[45]],
	#				results[[46]], results[[47]], results[[48]]
	#					orderBy = 1, ranksOf = 1, topNodes = length(allGO), numChar = 100)
	allRes <- GenTable(sampleGOdata, results[[1]], results[[2]], 
	orderBy = 1, ranksOf = 1, topNodes = length(allGO), numChar = 100)
	colnames(allRes)[6:ncol(allRes)] <- names(results)
	for(h in 6:ncol(allRes)) allRes[allRes[,h] == "< 1e-30",h] = 1e-30 
	for(h in 6:ncol(allRes)) class(allRes[,h]) = "numeric" 
	for(i in 1:sum(algoComp)){
		goAnalysis <- names(results[i])
		outputPie <- paste("topGO_", goTerm,"_", goAnalysis, "_pieChart.pdf", sep = "")
		outputNodes <- paste("topGO_", goTerm,"_", goAnalysis, "_nodes.pdf", sep = "")
		goDfr <- allRes
		colnames(goDfr)[1] <- "GO"
		goDfr <- goDfr[goDfr[,which(colnames(allRes) == goAnalysis)] < 0.05,]
		writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder,  resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("topGO_",goTerm,"_", goAnalysis))		
		drawPie(allRes, outputName = outputPie, output_folder = output_folder, 
		        pvalCol = which(colnames(allRes) == goAnalysis), 
		        labCol = 2, countCol = 4, totalCol = 3)
		
		pdf(paste(output_folder,outputNodes, sep = "/"), width = 21, height = 21, pointsize = 5)
		showSigOfNodes(sampleGOdata, score(results[[i]]), firstSigNodes = 10, useInfo ='all')
		dev.off()
	}
	output = paste(output_folder,"/topGO_",goTerm,"_allResults.tab",sep = "")
	write.table(allRes, output,sep = "\t")
	#return(allRes)
}
