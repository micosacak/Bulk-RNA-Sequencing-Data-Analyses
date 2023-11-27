runClusterProfiler <- function(ensembUniverse = NULL, selGenes = NULL, dfMatrix = NULL, 
                               output_folder = NULL, orgInfo = NULL, degDfr = NULL, 
                               geneIDtype = "ENSEMBL", prefix = "clusterProfiler", pvalCutoff = 0.05, ...){
  require("clusterProfiler")
  orgDb = eval(parse(text = orgInfo$orgDb))
  for(goType in c("BP","CC","MF")){
    enrichGOres <- enrichGO(gene = selGenes, universe = ensembUniverse, keyType = "ENSEMBL",
                            OrgDb = orgDb, ont = goType, pAdjustMethod = "BH",
                            pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
    enrichGOres = as.data.frame(enrichGOres)

    if(nrow(enrichGOres) > 0){  
            colnames(enrichGOres)[c(1,2,5,9)] = c("GO:ID","Term","Pvalue","Count")
    head(enrichGOres)
    enrichGOres$GeneRatio = ldply(sapply(enrichGOres$GeneRatio, strsplit, split = "/"))[,2]
    enrichGOres$BgRatio = ldply(sapply(enrichGOres$BgRatio, strsplit, split = "/"))[,2]

    drawPie(Df = enrichGOres,
            pvalCutoff = pvalCutoff,
            outputName = paste("CP_enrichGO_",goType,".pdf", sep = ""), output_folder = output_folder,
            pvalCol = 5,
            labCol = 2,
            countCol = 9,
            totalCol = 4)
    writeGOs(goDfr = enrichGOres, orgDb = orgInfo$orgDb, output_folder = output_folder,
             resultDfr = degDfr, geneIDtype = geneIDtype, goType = paste0("CP_",goType))
    }

  }
  
  entGenes = select(orgDb, keys = selGenes, columns = c("ENTREZID"), keytype = c("ENSEMBL"))
  entGenes = entGenes$ENTREZID
  entGenes = entGenes[!is.na(entGenes)]
  keggRes <- enrichKEGG(gene = entGenes, organism = orgInfo$org3L, pvalueCutoff = 0.05)
  keggRes = as.data.frame(keggRes)  
  
  if(nrow(keggRes) > 0) {
  
  colnames(keggRes)[c(1,2,5,9)] = c("KEGGID","Term","Pvalue","Count")
  keggRes$KEGGID = gsub(paste("\\",orgInfo$org3L, sep = ""),"",keggRes$KEGGID)
  
  keggRes$GeneRatio = ldply(sapply(keggRes$GeneRatio, strsplit, split = "/"))[,2]
  keggRes$BgRatio = ldply(sapply(keggRes$BgRatio, strsplit, split = "/"))[,2]

  drawPie(Df = keggRes,
          pvalCutoff = pvalCutoff, 
          outputName = paste("CP_enrichKEGGS",".pdf", sep = ""), output_folder = output_folder, 
          pvalCol = 5, 
          labCol = 2, 
          countCol = 9, 
          totalCol = 4)	
  
  drawKEGG(dfMatrix = dfMatrix, keggids = keggRes, output_folder = output_folder, gene.idtype = geneIDtype, keggNative = TRUE)
  
  writePathways(pathIds = gsub(orgInfo$org3L,"",keggRes$KEGGID), orgDb = orgDb, output_folder = output_folder, 
                resultDfr = degDfr, geneIDtype = geneIDtype, prefix = prefix)
  }
}
