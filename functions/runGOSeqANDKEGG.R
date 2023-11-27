#updated on 20.10.2017
runGOSeqANDKEGG <- function(do_GO_analysis = TRUE, uniGenes = NULL, orgId = "Dr", output_folder = NULL, dfMatrix = NULL, samplings = 10000, pvalCutoff = 0.1,  degDfr = NULL, geneIDtype = "ENSEMBL", ...){
	#dir.create(paste(output_folder,"GOSeq", sep = "/") ,showWarnings = FALSE)
	if(!(orgId %in% c("Dr","Mm","Hs","Rn"))){
		stop("orgId must be 'Dr', 'Mm' or 'Hs'")
	}
	sqlite_file <- paste0(orgInfo$org3L,"_txdb.sqlite")
	orgDb <- eval(parse(text = orgInfo$orgDb))
	if(sqlite_file %in% dir("rdaFiles/")){
		print(paste(sqlite_file," exists! Loading from the file!"), sep = "")
		txdb <- loadDb(paste("rdaFiles/",sqlite_file, sep = ""))
	}else{
		print(paste(sqlite_file," does not exist! Loading from the BiomaRt!"), sep = "")
		txdb <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = orgInfo$orgEs, transcript_ids = NULL, circ_seqs = DEFAULT_CIRC_SEQS, id_prefix = "ensembl_", 
		host = "www.ensembl.org", port = 80, taxonomyId = NA, miRBaseBuild = NA)
		saveDb(txdb,paste("rdaFiles/",sqlite_file, sep = ""))
	}
	txsByGene <- transcriptsBy(txdb,"gene")
	lengthData <- median(width(txsByGene))	
	uniGenes <- uniGenes[names(uniGenes) %in% names(lengthData)]
	lengthData <- lengthData[names(lengthData) %in% names(uniGenes)]	
	print("Now calculating KEGG Pathway")
	ensID2pathID <- as.list(sapply((mapIds(orgDb,  keys = keys(orgDb, keytype = "ENSEMBL"),column = "PATH",  keytype = "ENSEMBL", multiVals = "CharacterList")), paste, collapse = ","))
	pwf <- nullp(uniGenes, orgInfo$orgGb, "ensGene", bias.data = lengthData)
	GO.kegg <- goseq(pwf, gene2cat = ensID2pathID, use_genes_without_cat = TRUE )
	write.table(GO.kegg,paste(output_folder,paste("GO.kegg.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
	adjusted.GO.kegg.Over.BH <- GO.kegg[p.adjust(GO.kegg$over_represented_pvalue,method="BH") < pvalCutoff,]	
	adjusted.GO.kegg.Over.BH <- adjusted.GO.kegg.Over.BH[!(adjusted.GO.kegg.Over.BH$category == "NA"),]
	if(nrow(adjusted.GO.kegg.Over.BH) != 0){
		adjusted.GO.kegg.Over.BH <- cSplit(adjusted.GO.kegg.Over.BH, "category", ",", "long") # splits the category values, and generates a new row!
		adjusted.GO.kegg.Over.BH$category <- str_pad(adjusted.GO.kegg.Over.BH$category , 5, pad = "0")
		adjusted.GO.kegg.Over.BH <- adjusted.GO.kegg.Over.BH[order(adjusted.GO.kegg.Over.BH$category),]
		adjusted.GO.kegg.Over.BH <- adjusted.GO.kegg.Over.BH[!duplicated(adjusted.GO.kegg.Over.BH$category),] # removes duplicated pathids!
		write.table(adjusted.GO.kegg.Over.BH, paste(output_folder,paste("adjusted.GO.kegg.Over.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
		pathNames <- as.data.frame(KEGGPATHID2NAME)
		pathNames <- pathNames[pathNames$path_id %in% adjusted.GO.kegg.Over.BH$category,]
		adjusted.GO.kegg.Over.BH <- merge(adjusted.GO.kegg.Over.BH, pathNames, by.x = "category", by.y = "path_id")
		colnames(adjusted.GO.kegg.Over.BH) <- c("KEGGID","Pvalue", "under_represented_pvalue","Count", "numInCat","Term")
		if(startsWith(rownames(dfMatrix)[1], "ENS")){
			gene.idtype = "ENSEMBL"
		}else{
			gene.idtype = "ENTREZ"
		}
		drawPie(Df = as.data.frame(adjusted.GO.kegg.Over.BH), pvalCutoff = pvalCutoff, outputName = paste("KEGG","adjusted.GO.kegg.Over.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 2,labCol = 6, countCol = 4, totalCol = 5)
		if(nrow(dfMatrix) != 0 & ncol(dfMatrix) != 0){
		 	drawKEGG(dfMatrix, keggids = adjusted.GO.kegg.Over.BH, output_folder = output_folder, gene.idtype)
		}	
		pathIds <- unique(adjusted.GO.kegg.Over.BH$KEGGID)
		print(pathIds)
		writePathways(pathIds = pathIds, orgDb = eval(parse(text = paste("org",orgId,"eg.db", sep = "."))), output_folder, resultDfr = degDfr, geneIDtype = geneIDtype, prefix = "GOSeq")
	}	
	# adjusted.GO.kegg.Under.BH <- GO.kegg[p.adjust(GO.kegg$under_represented_pvalue,method="BH") < pvalCutoff,]
	# adjusted.GO.kegg.Under.BH <- adjusted.GO.kegg.Under.BH[!(adjusted.GO.kegg.Under.BH$category == "NA"),]
	# if(nrow(adjusted.GO.kegg.Under.BH) != 0){
	# 	adjusted.GO.kegg.Under.BH <- cSplit(adjusted.GO.kegg.Under.BH, "category", ",", "long")
	# 	adjusted.GO.kegg.Under.BH$category <- str_pad(adjusted.GO.kegg.Under.BH$category , 5, pad = "0")
	# 	write.table(adjusted.GO.kegg.Under.BH, paste(output_folder,paste("adjusted.GO.kegg.Under.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
	# 	pathNames <- as.data.frame(KEGGPATHID2NAME)
	# 	pathNames <- pathNames[pathNames$path_id %in% adjusted.GO.kegg.Under.BH$category,]
	# 	adjusted.GO.kegg.Under.BH <- merge(adjusted.GO.kegg.Under.BH, pathNames, by.x = "category", by.y = "path_id")
	# 	colnames(adjusted.GO.kegg.Under.BH) <- c("KEGGID","Pvalue", "under_represented_pvalue","Count", "numInCat","Term")
	# 	if(startsWith(rownames(dfMatrix)[1], "ENS")){
	# 		gene.idtype = "ENSEMBL"
	# 	}else{
	# 		gene.idtype = "ENTREZ"
	# 	}
	# 	drawPie(Df = as.data.frame(adjusted.GO.kegg.Under.BH), pvalCutoff = pvalCutoff, outputName = paste("KEGG","adjusted.GO.kegg.Under.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 2,labCol = 6, countCol = 4, totalCol = 5)
	# 	if(nrow(dfMatrix) != 0 & ncol(dfMatrix) != 0){
	# 		drawKEGG(dfMatrix, keggids = adjusted.GO.kegg.Under.BH, output_folder = output_folder, gene.idtype)
	# 	}	
	# 	pathIds <- unique(adjusted.GO.kegg.Under.BH$KEGGID)
	# 	writePathways(pathIds = pathIds, orgDb = eval(parse(text = paste("org",orgId,"eg.db", sep = "."))), output_folder, resultDfr = degDfr, geneIDtype = geneIDtype, prefix = "GOSeq")
	# }	
	# 
	pwf = nullp(uniGenes,orgInfo$orgGb,"ensGene", bias.data = lengthData)
	# choose over- and under-represented categories with adjusting p-values!
	#p.adjust(p, method = p.adjust.methods, n = length(p))
	#p.adjust.methods // c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")//	
	if(do_GO_analysis){
		for(goType in c("BP","CC","MF")){
			print(goType)
			GO.wall = goseq(pwf,orgInfo$orgGb,"ensGene", method = "Wallenius", use_genes_without_cat = TRUE,test.cats = c(paste("GO",goType, sep = ":")))
			GO.samp = goseq(pwf,orgInfo$orgGb,"ensGene", method = "Sampling",repcnt=samplings, use_genes_without_cat = TRUE, test.cats = c(paste("GO",goType, sep = ":")))
			GO.nobias = goseq(pwf,orgInfo$orgGb,"ensGene",method = "Hypergeometric",use_genes_without_cat = TRUE, test.cats = c(paste("GO",goType, sep = ":")))
			pdf(paste(output_folder,paste(goType,"Wallenius_vs_Sampling.pdf", sep = "_"), sep ="/"),height = 14, width =14, pointsize = 12)
			plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]),
			xlab="log10(Wallenius p-values)",ylab = "log10(Sampling p-values)",xlim = c(-3,0))
			abline(0,1,col=3,lty=2)
			dev.off()
			pdf(paste(output_folder,paste(goType,"Wallenius_vs_Hypergeometic.pdf", sep = "_"), sep ="/"),height = 14, width = 14, pointsize = 12)
			plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.nobias[,1],GO.wall[,1]),2]),
			xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)",xlim=c(-3,0), ylim=c(-3,0))
			abline(0,1,col=3,lty=2)
			dev.off()
			pdf(paste(output_folder,paste(goType,"Sampling_vs_Hypergeometic.pdf", sep = "_"), sep ="/"),height = 14, width = 14, pointsize = 12)
			plot(log10(GO.samp[,2]), log10(GO.nobias[match(GO.nobias[,1],GO.samp[,1]),2]),
			xlab="log10(Sampling p-values)", ylab="log10(Hypergeometric p-values)",xlim=c(-3,0), ylim=c(-3,0))
			abline(0,1,col=3,lty=2)
			dev.off()
			adjusted.GO.wall.Over.BH <- GO.wall[p.adjust(GO.wall$over_represented_pvalue,method="BH") < pvalCutoff,]
			adjusted.GO.wall.Under.BH <- GO.wall[p.adjust(GO.wall$under_represented_pvalue,method="BH") < pvalCutoff,]
			adjusted.GO.samp.Over.BH <- GO.samp[p.adjust(GO.samp$over_represented_pvalue,method="BH") < pvalCutoff,]
			adjusted.GO.samp.Under.BH <- GO.samp[p.adjust(GO.samp$under_represented_pvalue,method="BH") < pvalCutoff,]	
			adjusted.GO.nobias.Over.BH <- GO.nobias[p.adjust(GO.nobias$over_represented_pvalue,method="BH") < pvalCutoff,]	
			adjusted.GO.nobias.Under.BH <- GO.nobias[p.adjust(GO.nobias$under_represented_pvalue,method="BH") < pvalCutoff,]
			write.table(GO.wall,paste(output_folder,paste(goType,"GO.wall.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(GO.samp,paste(output_folder,paste(goType,"GO.samp.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(GO.nobias,paste(output_folder,paste(goType,"GO.nobias.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(adjusted.GO.wall.Over.BH,paste(output_folder,paste(goType,"adj.GO.wall.Over.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(adjusted.GO.wall.Under.BH,paste(output_folder,paste(goType,"adj.GO.wall.Under.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(adjusted.GO.samp.Over.BH,paste(output_folder,paste(goType,"adj.GO.samp.Over.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(adjusted.GO.samp.Under.BH,paste(output_folder,paste(goType,"adj.GO.samp.Under.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(adjusted.GO.nobias.Over.BH,paste(output_folder,paste(goType,"adj.GO.nobias.Over.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")
			write.table(adjusted.GO.nobias.Under.BH,paste(output_folder,paste(goType,"adj.GO.nobias.Under.BH.tab", sep = "_"), sep = "/"), row.names = FALSE, sep = "\t")

			goDfr <- adjusted.GO.wall.Over.BH
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOSeq_",goType,".adj.GO.wall.Over.BH"))
			goDfr <- adjusted.GO.wall.Under.BH
			colnames(goDfr)[c(1,6)] <- c("GO","Term")
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOSeq_",goType,".adj.GO.wall.Under.BH"))
			goDfr <- adjusted.GO.samp.Over.BH
			colnames(goDfr)[c(1,6)] <- c("GO","Term")
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOSeq_",goType,".adj.GO.samp.Over.BH"))
			goDfr <- adjusted.GO.samp.Under.BH
			colnames(goDfr)[c(1,6)] <- c("GO","Term")
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOSeq_",goType,".adj.GO.samp.Under.BH"))
			goDfr <- adjusted.GO.nobias.Over.BH
			colnames(goDfr)[c(1,6)] <- c("GO","Term")
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOSeq_",goType,".adj.GO.nobias.Over.BH"))
			goDfr <- adjusted.GO.nobias.Under.BH
			colnames(goDfr)[c(1,6)] <- c("GO","Term")
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOSeq_",goType,".adj.GO.nobias.Under.BH"))

			drawPie(Df = adjusted.GO.wall.Over.BH, pvalCutoff = 0.1, outputName = paste(goType,"adj.GO.wall.Over.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 2,labCol = 6 , countCol = 4, totalCol = 5)
			drawPie(Df = adjusted.GO.wall.Under.BH, pvalCutoff = 0.1, outputName = paste(goType,"adj.GO.wall.Under.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 3,labCol = 6 , countCol = 4, totalCol = 5)
			drawPie(Df = adjusted.GO.samp.Over.BH, pvalCutoff = 0.1, outputName = paste(goType,"adj.GO.samp.Over.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 2,labCol = 6 , countCol = 4, totalCol = 5)
			drawPie(Df = adjusted.GO.samp.Under.BH, pvalCutoff = 0.1, outputName = paste(goType,"adj.GO.samp.Under.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 3,labCol = 6 , countCol = 4, totalCol = 5)
			drawPie(Df = adjusted.GO.nobias.Over.BH, pvalCutoff = 0.1, outputName = paste(goType,"adj.GO.nobias.Over.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 2,labCol = 6 , countCol = 4, totalCol = 5)
			drawPie(Df = adjusted.GO.nobias.Under.BH, pvalCutoff = 0.1, outputName = paste(goType,"adj.GO.nobias.Under.BH.pdf", sep = "_"), output_folder = output_folder, pvalCol = 3,labCol = 6 , countCol = 4, totalCol = 5)
		}	
	}
		# save the results!
}
