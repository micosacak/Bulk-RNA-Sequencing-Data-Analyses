#updated on 20.10.2017
runGOStats <- function(do_GO_analysis = FALSE, orgId = "Dr", uniGenes = NULL, selGenes = NULL, output_folder = NULL, dfMatrix = NULL, pvalCutoff = 0.05, geneIDtype = "ENTREZ", organism = "Danio rerio", degDfr = NULL,...){	
	dir.create(output_folder ,showWarnings = FALSE)
	files = dir("rdaFiles")
	if(!(orgId %in% c("Dr","Mm","Hs","Rn"))){
		stop("orgId must be 'Dr', 'Mm' or 'Hs'")
	}
	if(geneIDtype == "ENTREZ"){
		efileName <- paste("gsc",orgId,"entrezids.rda", sep = "")
		kfileName <- paste("kgsc",orgId,"entrezids.rda", sep = "")
	}else if(geneIDtype == "ENSEMBL"){
		efileName <- paste("gsc",orgId,"ensemblids.rda", sep = "")
		kfileName <- paste("kgsc",orgId,"ensemblids.rda", sep = "")
	}else{
		stop("gene id type must be ENSEMBL or ENTREZ")
	}
	egENS <- keys(eval(parse(text = paste("org",orgId,"eg.db", sep = "."))))
	egENS <- AnnotationDbi::select(eval(parse(text = paste("org",orgId,"eg.db", sep = "."))), keys = egENS, columns = "ENSEMBL", keytype = "ENTREZID")
	if(!(efileName %in% files)){
		dframe <- getGOframe(geneIDtype, eval(parse(text = paste("org",orgId,"egGO", sep = "."))), eval(parse(text = paste("org",orgId,"eg.db", sep = "."))))
		goframeData = data.frame(dframe$go_id, dframe$Evidence, dframe$gene_id)
		goFrame = GOFrame(goframeData, organism = organism)
		goAllFrame = GOAllFrame(goFrame)
		gscOrg <- GeneSetCollection(goAllFrame, setType = GOCollection())
		save(gscOrg, file = paste("rdaFiles/",efileName, sep = ""))
	}else{
		load(file = paste("rdaFiles/",efileName, sep = ""))
	}
	if(!(kfileName %in% files)){
		kframe = toTable(eval(parse(text = paste("org",orgId,"egPATH", sep = "."))))
		keggFrameData = data.frame(kframe$path_id, kframe$gene_id)
		keggFrame=KEGGFrame(keggFrameData,organism = organism)
		kgscOrg <- GeneSetCollection(keggFrame, setType = KEGGCollection())
		save(kgscOrg, file = paste("rdaFiles/",kfileName, sep = ""))
	}else{
		load(file = paste("rdaFiles/",kfileName, sep = ""))
	}
	goTerms <- c("MF","CC","BP")
	result <- list("MF_Up" = NA,"MF_Down" = NA,"CC_Up" = NA,"CC_Down" = NA,"BP_Up" = NA,"BP_Down" = NA, "kegg_Up" = NA, "kegg_Down" = NA)
	GOstatssm <- GOstats::summary
	if(do_GO_analysis){
		for(goTerm in goTerms){
			hypTest <- paste(goTerm,"_Up", sep = "")
			print(hypTest)
			params <- GSEAGOHyperGParams(name="Custom en2GO GSEABase",
				geneSetCollection = gscOrg,
				geneIds = selGenes,
				universeGeneIds = uniGenes,
				ontology = goTerm,
				pvalueCutoff = pvalCutoff,
				conditional = FALSE,
				testDirection = "over")
			Over <- hyperGTest(params)
			result[[which(names(result) == hypTest)]] = summary(Over)

			goDfr <- summary(Over)
			colnames(goDfr)[1] <- "GO"
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOstats_",hypTest))
			drawPie(Df = summary(Over), pvalCutoff = pvalCutoff, outputName = paste("GOstats_",hypTest,"_pieChart.pdf", sep = ""), output_folder = output_folder, pvalCol = 2, labCol = 7, countCol = 5, totalCol = 6)
			write.table(summary(Over), paste(output_folder, paste("GOstats_",hypTest,".tab", sep = ""), sep = "/"), sep = "\t")
			
			hypTest <- paste(goTerm,"_Down", sep = "")
			print(hypTest)
			params <- GSEAGOHyperGParams(name="Custom en2GO GSEABase",
				geneSetCollection = gscOrg,
				geneIds = selGenes,
				universeGeneIds = uniGenes,
				ontology = goTerm,
				pvalueCutoff = pvalCutoff,
				conditional = FALSE,
				testDirection = "under")
			Under <- hyperGTest(params)
			result[[which(names(result) == hypTest)]] = summary(Under)

			goDfr <- summary(Over)
			colnames(goDfr)[1] <- "GO"
			writeGOs(goDfr = goDfr, orgDb = orgInfo$orgDb, output_folder = output_folder, resultDfr = degDfr, geneIDtype = "ENSEMBL", goType = paste0("GOstats_",hypTest))

			drawPie(Df = summary(Under), pvalCutoff = pvalCutoff, outputName = paste("GOstats_",hypTest,"_pieChart.pdf", sep = ""), output_folder = output_folder, pvalCol = 2, labCol = 7, countCol = 5, totalCol = 6)
			write.table(summary(Under), paste(output_folder, paste("GOstats_",hypTest,".tab", sep = ""), sep = "/"), sep = "\t")
		}
	}
	
	# perform kegg pathway analysis.
	if(geneIDtype == "ENSEMBL"){
		uniGenes <- egENS$ENTREZID[egENS$ENSEMBL %in% uniGenes]
		selGenes <- egENS$ENTREZID[egENS$ENSEMBL %in% selGenes]
	}
	print("kegg_Up")
	kparams <- GSEAKEGGHyperGParams(name="Custom en2kegg GSEABase",
		geneSetCollection= kgscOrg,
		geneIds = selGenes,
		universeGeneIds = uniGenes,
		pvalueCutoff = pvalCutoff,
		testDirection = "over")
	kOver <- hyperGTest(kparams)	
	pathWays <- summary(kOver)
	pathIds <- pathWays$KEGGID
	writePathways(pathIds = pathIds, orgDb = eval(parse(text = paste("org",orgId,"eg.db", sep = "."))), output_folder = output_folder, resultDfr = degDfr, geneIDtype = geneIDtype, prefix = "GOstats")
	drawPie(Df = summary(kOver), pvalCutoff = pvalCutoff, outputName = paste("GOstats_","kegg_Up",".pdf", sep = ""), output_folder = output_folder, pvalCol = 2, labCol = 7, countCol = 5, totalCol = 6)	
	
	print(rownames(dfMatrix[1]))
	print(head(dfMatrix))
	if(startsWith(rownames(dfMatrix)[1], "ENS")){
		gene.idtype = "ENSEMBL"
	}else{
		gene.idtype = "ENTREZ"
	}
	if(nrow(dfMatrix) != 0 & ncol(dfMatrix) != 0){
		drawKEGG(dfMatrix, summary(kOver), output_folder, gene.idtype)
	}		
	result[[which(names(result) == "kegg_Up")]] = summary(kOver)
	write.table(summary(kOver), paste(output_folder, paste("GOstats_kegg_over.tab", sep = ""), sep = "/"), sep = "\t")
	print("kegg_Down")
	kparams <- GSEAKEGGHyperGParams(name="Custom en2kegg GSEABase",
		geneSetCollection= kgscOrg,
		geneIds = selGenes,
		universeGeneIds = uniGenes,
		pvalueCutoff = pvalCutoff,
		testDirection = "under")
	kUnder <- hyperGTest(kparams)	
	result[[which(names(result) == "kegg_Down")]] = summary(kUnder)
	drawPie(Df = summary(kUnder), pvalCutoff = pvalCutoff, outputName = paste("GOstats_","kegg_Under",".pdf", sep = ""), output_folder = output_folder, pvalCol = 2, labCol = 7, countCol = 5, totalCol = 6)		
	write.table(summary(kUnder), paste(output_folder, paste("GOstats_kegg_under.tab", sep = ""), sep = "/"), sep = "\t")
}
