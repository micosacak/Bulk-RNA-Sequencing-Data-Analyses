# 04.11.2017

# load all functions in fucntions that might be required!
functions = dir("functions/")
for(func in functions) source(paste("functions/", func, sep = ""))
functions = NULL
if(("heatmap3_function.R" %in% dir("functions/"))){
  source("functions/heatmap3_function.R")
}else{
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
}

# install and load libraries
installLoadLibs()
select = AnnotationDbi::select  # make sure that "select" is from "AnnotationDbi" not from "MASS"!!!

# some default values
savelogFCpadj = TRUE     # saves the number of DE genes based on the log2 threshold
saveAllData = TRUE       # saves rawData & normalization data

fdr_threshold = 0.1          # FDR or p.adjust Thresholds
lfc_threshold = log2(2)      # log2 Fold Change Threshold

mainDir <<- getwd()             # global required!!!
theme_set(theme_bw(16))         # set theme!

dir.create("rdaFiles", showWarnings = FALSE)

# set sysmtem info! check if it is Windows or Linux. register cores and also endswith function!
if(Sys.info()["sysname"][[1]] == "Linux"){
  endsWith <- stringi::stri_endswith # endsWith does not work on Linux!!! instead, use stri_endswith
  register(MulticoreParam(detectCores()-2))
}else if(Sys.info()["sysname"][[1]] == "Windows"){
  register(SnowParam(detectCores()-2))
}

print("Now getting the inputs ...")
#### The Input Files must have column names as Sample1_R1, Sample1_R2, ...
#### The default input file is from featureCounts!
#### other wise indicate which columns to be removed!
#### which column is rownames!
#### the rmSuf must be indicated correctly! (the suffic to be removed)
#### OR user must provide colnames seperately
if(!("inputs.rda" %in% dir("rdaFiles"))){
  input_files = dir("inputFiles")
  input_files = input_files[endsWith(input_files, ".tab")]
  inputRDAs = vector("list", length = length(input_files))
  for(i in 1:length(input_files)) inputRDAs[[i]] = readInputFiles(input_file = paste("inputFiles/",input_files[i], sep = ""), rmSuf = ".fastq.gz.sam")
  rawCounts = mergeAllDfr(inputRDAs)
  colnames(rawCounts) = gsub("L\\d+_","",colnames(rawCounts))        #### if required, change colnames
  colnames(rawCounts)
  all_conditions = gsub("\\_R[0-9]","",colnames(rawCounts))
  print(unique(all_conditions))
  comparisons = c(1,2, 3,4, 6,5, 4,5, 4,2) 
  print(comparisons)
  # check all_conditions and conditions to decide comparisons
  inputs = list("rawCounts" = rawCounts, "comps" = comparisons,"all_conditions" = all_conditions)
  save(inputs, file = "rdaFiles/inputs.rda")
}else{
  print("loading inpus.rda file")
  load("rdaFiles/inputs.rda")
}

# some more default values
orgId <- "Hs"
orgInfo <<- getOrgInfo(orgId) # orgInfo as global object!
ensembUniverse <- unique(rownames(inputs$rawCounts))
entrezUniverse <- unique(as.character(keys(eval(parse(text = orgInfo$orgDb)))))

all_conditions = inputs$all_conditions
comps = unique(all_conditions)
aCo = c()
for(i in seq(2,length(inputs$comps),2)){
  aCo = c(paste0(comps[inputs$comps[i-1]],"_vs_",comps[inputs$comps[i]]),aCo)
}
print(aCo)   # check comparisons!

outFile = "results_on"
tm = as.character(Sys.time())
tm <- gsub("[-: ]","",tm)
outFile1 <- paste0(outFile,"_",tm)
dir.create(outFile1, showWarnings = FALSE)

### start analysis
print("Now Starting Analysis ...")

if(!(paste0(outFile,".rda") %in% dir("rdaFiles/")) || !("dds.rda" %in% dir("rdaFiles/"))){
  print("Now generating rda file ...")
  coldata <- data.frame(row.names = colnames(inputs$rawCounts), condition = as.factor(all_conditions))
  dds <- DESeqDataSetFromMatrix(countData = inputs$rawCounts, colData = coldata, design = ~ condition)
  dds <- estimateSizeFactors(dds)
  dds = DESeq(dds, parallel = T)
  rsDfr <- getDESeq2Results(aCo = aCo, dds = dds, lfc_threshold = lfc_threshold, fdr_threshold = fdr_threshold) # added on 22.09.2016!
  save(rsDfr,file = paste0("rdaFiles/",outFile,".rda"))
  save(dds, file = "rdaFiles/dds.rda")
}else{
  print("Now loading rda file ...")
  load("rdaFiles/dds.rda")
  load(paste0("rdaFiles/",outFile,".rda", sep = ""))
}	

# as user can choose different threshold for significantly expressed genes
# it is better to show genes with three options:
# 0  : no change in the expression level; thus abs(log2FoldChange) < lfc_threshold
# -1 : down-regulated genes; thus log2FoldChange <= lfc_threshold & padj < fdr_threshol
# 1 : up-regulated genes; thus log2FoldChange >= lfc_threshold & padj < fdr_threshol
# the heatmap below is just using 0: no change, -1: down-regulation and +1:upregulation
theMatrix <- rsDfr$dfr$deseq$pathDfr
theMatrix <- theMatrix[rowSums(abs(theMatrix)) != 0,]
#drawHeatMap(theMatrix = as.matrix(theMatrix), output_folder = paste0(mainDir,"/",outFile1))

# plot the number of differentially expressed genes versus log2FoldChange!
# 
# if(savelogFCpadj){
#   deseq_cols <- sort(c(seq(2,ncol(rsDfr$dfr$deseq$outDfr),6),seq(6,ncol(rsDfr$dfr$deseq$outDfr),6)))
#   rsDfr$dfr$deseq$outDfr <- rsDfr$dfr$deseq$outDfr[,deseq_cols]
#   plotLogs(orgDfr = rsDfr$dfr$deseq$outDfr, output_folder = paste0(mainDir,"/",outFile1,"/","deseq2_logFCs.pdf", sep = ""))
# }	

# better to save raw counts, normalized counts and DESeq2 results

if(saveAllData){
  rsDfr$dfr$deseq$outDfr <- mergeAllDfr(list("1" = rsDfr$dfr$deseq$outDfr, "2" = inputs$rawCounts, "3" = as.data.frame(counts(dds, normalized = TRUE))))
}

outFile = outFile1
nComp = 1
laCo <- length(aCo)
showAll = FALSE
i = 1
rs <- rsDfr$rs
dfr <- rsDfr$dfr

drawVenn(dds = dds, aCo = aCo, output_folder = outFile)

drawVenn(dds = dds, aCo = aCo[3:5], output_folder = outFile)

if(!("rld.rda" %in% dir(mainDir))){
  rld <<- rlogTransformation(dds)
  save(rld, file =  "rdaFiles/rld.rda")
}else{
  load("rdaFiles/rld.rda")
}
plotSampleDists(rld, outFile, all_conditions)


pdf(paste(outFile,"PCA_allsamples.pdf", sep = "/"))
print(plotPCA(rld, intgroup = "condition"))
dev.off()

# get annotaions: SYMBOL and GENENAME
orgDb = eval(parse(text = paste("org",orgId,"eg.db", sep = ".")))
ann = select(orgDb, keys = rownames(inputs$rawCounts), columns = c("SYMBOL","GENENAME"), keytype = "ENSEMBL")
ann = ann[!duplicated(ann$ENSEMBL),]
rownames(ann) = ann$ENSEMBL
ann$ENSEMBL = NULL
head(ann)
normCounts = counts(dds, normalized = TRUE)
rawCounts = counts(dds, normalized = FALSE)
allResults = mergeAllDfr(list(rawCounts,normCounts,ann))


for(i in 1:length(aCo)){
  output_folder <- aCo[i]
  outputFolder <- paste0(outFile,"/",output_folder, collapse = "")
  dir.create(outputFolder ,showWarnings = FALSE)
  
  print(paste("Now Analysing  ",output_folder," ...", sep = ""))
  dir.create(paste(outputFolder, sep = "/"), showWarnings = FALSE)
  output_file_name <- paste(output_folder, "diffexp-results.xlsx", sep = "_") 
  
  dds_res <- rs$deseq[[i]]
  resdata <- mergeAllDfr(list(dds_res,allResults))
  chooseCols = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj",unlist(strsplit(output_folder,"_vs_")),"SYMBOL","GENENAME")
  outDfr = resdata[,c(gsub("_R[0-9].[x-y]","",colnames(resdata)) %in% chooseCols)]
  write.xlsx(outDfr, sheetName = output_folder, file = paste(outputFolder,paste("allResults",output_file_name,sep = "_"), sep = "/"), row.names = T, col.names = T)
  sig_resdata_deg <- subset(outDfr, outDfr$padj <= fdr_threshold & abs(outDfr$log2FoldChange) >= lfc_threshold)
  write.xlsx(sig_resdata_deg,sheetName = output_folder, file = paste(outputFolder,paste("sigDEGs",output_file_name,sep = "_"), sep = "/"), row.names = T, col.names = T)
  write.xlsx(list("allResults"= outDfr,"sigDEGs" = sig_resdata_deg), file = paste(outputFolder,paste("all_and_sig",output_file_name,sep = "_"), sep = "/"), row.names = T, col.names = T)
  
  
  pdf(paste(outputFolder,paste(output_folder,"_MAplot.pdf", sep = ""), sep = "/"))
  plotMASmears(resdata[,c(1,2,6)], main = output_folder, LoG = TRUE, lfc_threshold = lfc_threshold, xNames = "log10basemean")
  dev.off()
  
  pdf(paste(outputFolder,paste(output_folder,"_Volcanoplot.pdf", sep = ""), sep = "/"))
  plotVolcano(resdata[,c(2,5,6)], main = output_folder, lfc_threshold = lfc_threshold)
  dev.off()
  
  dfm <- getDfMatrix(resdata[,c(2,6)], lfc_threshold = lfc_threshold)
  
  runGOpath = TRUE
  universe <- ensembUniverse
  selGenes <- unique(rownames(dfm)[abs(dfm$logFC) == 1])
  uniGenes <- as.integer(universe %in% selGenes)
  names(uniGenes) <- universe
  
  if(showAll){
    dfMatrixDfr = data.frame(dfr$deseq$pathDfr[,nComp:laCo], row.names = rownames(dfr$deseq$pathDfr))
    colnames(dfMatrixDfr) = aCo[nComp:laCo]
  }else{
    dfMatrixDfr = data.frame(dfr$deseq$pathDfr[,i:i], row.names = rownames(dfr$deseq$pathDfr))
    colnames(dfMatrixDfr) = aCo[i:i]
  }
  
  if(runGOpath){
    dir.create(paste(outputFolder,"clusterProfiler", sep = "/"), showWarnings = FALSE)
    runClusterProfiler(ensembUniverse = ensembUniverse, selGenes = selGenes, dfMatrix = dfMatrixDfr, 
                       output_folder = paste(outputFolder,"clusterProfiler", sep = "/"), 
                       orgInfo = orgInfo, degDfr = outDfr)
    
    dir.create(paste(outputFolder,"GOstats", sep = "/"), showWarnings = FALSE)
    runGOStats(orgId = orgId, uniGenes = universe, selGenes = selGenes, geneIDtype = "ENSEMBL", output_folder = paste(outputFolder,"GOstats", sep = "/"), 
               dfMatrix = dfMatrixDfr, 
               pvalCutoff = 0.05, organism = orgInfo$orgSp, degDfr = outDfr, do_GO_analysis = TRUE) # !!!!
    
    dir.create(paste(outputFolder,"GOSeq", sep = "/"), showWarnings = FALSE)
    runGOSeqANDKEGG(uniGenes = uniGenes, orgId = orgInfo$orgId, 
                    output_folder = paste(outputFolder,"GOSeq", sep = "/"), samplings = 10000, 
                    dfMatrix = dfMatrixDfr, do_GO_analysis = FALSE,
                    degDfr  = outDfr, geneIDtype = "ENSEMBL")
    
    
    dir.create(paste(outputFolder,"deseq2","topGO", sep = "/"), showWarnings = FALSE)
    topGOanalysis(uniGenes = uniGenes, output_folder = paste(outputFolder,"topGO", sep = "/"), degDfr = outDfr)
  }	
}

runHeatMaps = TRUE
if(runHeatMaps){
  ### heatmaps
  clstrFolder = paste(outFile1,"clusters", sep = "/")
  dir.create(clstrFolder, showWarnings = FALSE)
  if(!("ann.rda" %in% dir("rdaFiles/"))){
    library("org.Hs.eg.db")
    ann = select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = c("SYMBOL","GENENAME","ENSEMBL"))
    ann = ann[!is.na(ann$ENSEMBL),]
    ann = ann[!duplicated(ann$ENSEMBL),]
    rownames(ann) = ann$ENSEMBL
    ann = ann[,-c(1,3)]
    head(ann)
    save(ann, file = "rdaFiles/ann.rda")
  }else{
    load("rdaFiles/ann.rda")
  }
  
  if(!("candidateGenes.rda" %in% dir("rdaFiles/"))){
    candidateGenes = read.delim("inputFiles/candidateGenes.csv", header = T)
    candidateGenes = merge(candidateGenes,ann, by.x = "SYMBOL", by.y = "SYMBOL", all.x = TRUE, all.y = FALSE)
    candidateGenes = candidateGenes[!is.na(candidateGenes$ENSEMBL),]
    candidateGenes = candidateGenes[!duplicated(candidateGenes$ENSEMBL),]
    save(candidateGenes, file = "rdaFiles/candidateGenes.rda")
  }else{
    load("rdaFiles/candidateGenes.rda")
  }
  
  head(ann)
  head(candidateGenes)
  
  if(!("allNormData.rda" %in% dir("rdaFiles/"))){
    allNormData = allResults[,endsWith(colnames(allResults),"y")]
    colnames(allNormData) = gsub("\\.y","",colnames(allNormData))
    colnames(allNormData)
    df.mean = as.data.frame(matrix(0,nrow(allNormData),length(unique(all_conditions))))
    colnames(df.mean) = unique(all_conditions)
    rownames(df.mean) = rownames(allNormData)
    for(cN in colnames(df.mean)){
      idx = which(gsub("\\_R[0-9]","",colnames(allNormData)) == cN)
      df.mean[,cN] = rowMeans(allNormData[,idx])
    }
    df.mean$SYMBOL = q1$SYMBOL
    head(df.mean)
    allNormData = df.mean
    df.mean = NULL
    save(allNormData, file = "rdaFiles/allNormData.rda")
  }else{
    load("rdaFiles/allNormData.rda")
  }
  
  head(allNormData)
  
  if(!("mtrx.rda" %in% dir("rdaFiles/"))){
    mtrx = allNormData
    mtrx = mtrx[mtrx$SYMBOL %in% candidateGenes$SYMBOL,]
    mtrx = mtrx[!is.na(mtrx[,1]),]
    mtrx = merge(mtrx,candidateGenes,by = "SYMBOL", all = TRUE)
    mtrx = mtrx[!duplicated(mtrx$SYMBOL),]
    rownames(mtrx) = mtrx$SYMBOL
    mtrx$SYMBOL = NULL
    mtrx$ENSEMBL = NULL
    save(mtrx, file = "rdaFiles/mtrx.rda")
  }else{
    load("rdaFiles/mtrx.rda")
  }
  
  head(mtrx)
  
  mtrx = mtrx[mtrx$LAYER != "stemcell",]
  
  namesOrdered = c("RELN","MARCKSL1","PLXND1","SATB2",
                   "MDGA1","POU3F2","TLE3","LHX2","POU3F3",
                   "TLE1","CUX1","MEF2C","DTX4","SLITRK1",
                   "UNC5D","CUX2","RORB","KCNIP2","CYP39A1","S100A10",
                   "FOXG1","LIX1L","OMA1","CNTN6","RAC3","FOXO1","OPN3","LDB2","CRYM","DKK3","TLE4","SEMA3E",
                   "OTX1","SOX5","LXN","FOXP2","CRIM1","PCP4","CTGF","DLX1","IGFBP4","ID2")
  
  mtrx = mtrx[rownames(mtrx) %in% namesOrdered,]
  library("gplots")
  library("devtools")
  
  hclustfunc <- function(x) hclust(x, method="complete")
  distfunc <- function(x) dist(x,method="maximum")
  
  tmpMtrx = mtrx # keep for later use!
  
  tmpMat = mtrx
  tmpMat = tmpMat[rownames(tmpMat) %in% namesOrdered,]
  tmpMat = tmpMat[order(match(rownames(tmpMat),namesOrdered)),]
  
  mtrx  = tmpMat
  head(mtrx)
  mt = as.matrix(mtrx[,1:(ncol(mtrx)-1)])
  head(mt)
  
  getColors <- function(mtrx = NULL, colNames = c("1","2","3","4","5","6","subplate","generalcortex"), coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black")){
    myCols = as.data.frame(matrix(0,nrow(mtrx),length(colNames)))
    rownames(myCols) = rownames(mtrx)
    colnames(myCols) = colNames
    for(i in 1:nrow(mtrx)){
      ly = colNames
      sp = unlist(strsplit(as.character(mtrx$LAYER[i]),","))
      ly[!(ly %in% sp)] = 0
      myCols[i,] = ly
    }
    print(head(myCols))
    i = 6
    j = 3
    for(i in 1:nrow(myCols)){
      for(j in 1:ncol(myCols)){
        if(myCols[i,j] == 0){
          myCols[i,j] = "white"
        }else{
          myCols[i,j] = coloRs[which(colNames == myCols[i,j])]
        }
        
      }
    }
    colnames(myCols) = c("I","II","III","IV","V","VI","SP","G")
    return(as.matrix(myCols))
  }
  
  ### generate colors
  myCols = getColors(mtrx, colNames = c("1","2","3","4","5","6","subplate","generalcortex"), coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black"))
  
  head(myCols)
  mtrx$LAYER
  nrow(mtrx)
  nrow(myCols)
  
  tmpMt = mt
  
  getBreaksColors = function(maxLim = 500, lengthOut = 50000, minLim = 0){
    breaks = seq(minLim,max(mt),length.out= lengthOut)
    g1 = colorpanel(sum(breaks[-1] <= maxLim), "white","blue")
    g4 = colorpanel(sum(breaks[-1] > maxLim), "darkblue", "black")
    hm = c(g1,g4)
    return(hm)
  }
  getBreaksForKey = function(maxLim = 1000, rangeLim = 1001, lengthOut = 50000, minLim = 0){
    breaks = seq(minLim,max(maxLim),length.out= lengthOut)
    breaks[breaks >= rangeLim] = rangeLim
    return(breaks)
  }
  
  rangeLim = 1000
  colnames(mt)
  
  pdf(paste(clstrFolder,"/NoChangeInOrder_heatMap_allsamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4,RowSideColors = t(myCols), col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "Normalized Reads",Rowv = NA, Colv = NA, dendrogram = "none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(8,5))
  dev.off()
  
  colnames(mt)
  mt = mt[,-c(1,3,6)]
  colnames(mt)
  
  pdf(paste(clstrFolder,"/NoChangeInOrder_heatMap_ControlSamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4,RowSideColors = t(myCols), col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "Normalized Reads",Rowv = NA, Colv = NA, dendrogram = "none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(8,5))
  
  dev.off()
  
  # cluster namessss!
  getNameOrders = function(mtrx,colNames = c("1","2","3","4","5","6","subplate","generalcortex"),coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black")){
    myCols = getColors(mtrx, colNames = colNames, coloRs = coloRs)
    lv = levels(as.factor(myCols))
    myColMatrix = as.data.frame(matrix(0, nrow(myCols), ncol(myCols)))
    colnames(myColMatrix) = colnames(myCols)
    rownames(myColMatrix) = rownames(myCols)
    for(i in 1:ncol(myColMatrix)){
      for(j in 1:nrow(myColMatrix)){
        myColMatrix[j,i] = which(lv == myCols[j,i])
      }
    }
    head(myColMatrix)
    myColMatrix = as.matrix(myColMatrix)
    hh = heatmap.3(myColMatrix, hclustfun = hclustfunc, distfun = distfunc, RowSideColorsSize = 0.4,RowSideColors = t(myCols), key = F,
                   keysize = 0.9, KeyValueName = "Normalized Reads", dendrogram = "none",
                   cexRow = 0.7, cexCol = 1, lwid = c(0.2,0.2), margins = c(8,5))
    dev.off()
    hh$rowInd  
    return(rownames(myColMatrix)[hh$rowInd])
  }
  
  mtrx = tmpMtrx
  namesOrdered = getNameOrders(mtrx,colNames = c("1","2","3","4","5","6","subplate","generalcortex"),coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black"))
  tmpMat = mtrx
  tmpMat = tmpMat[rownames(tmpMat) %in% namesOrdered,]
  tmpMat = tmpMat[order(match(rownames(tmpMat),namesOrdered)),]
  
  mtrx  = tmpMat
  mt = as.matrix(mtrx[,1:(ncol(mtrx)-1)])
  head(mt)
  
  ### generate colors
  myCols = getColors(mtrx, colNames = c("1","2","3","4","5","6","subplate","generalcortex"), coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black"))
  
  head(myCols)
  mtrx$LAYER
  nrow(mtrx)
  nrow(myCols)
  
  tmpMt = mt
  
  colnames(mt)
  
  pdf(paste(clstrFolder,"/LayerClustered_heatMap_allsamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim),RowSideColorsSize = 0.4,RowSideColors = t(myCols), col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 0.9, KeyValueName = "Normalized Reads",Rowv = NA, Colv = NA, dendrogram = "none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(8,5))
  dev.off()
  
  colnames(mt)
  mt = mt[,-c(1,3,6)]
  colnames(mt)
  
  pdf(paste(clstrFolder,"/LayerClustered_heatMap_ControlSamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim),RowSideColorsSize = 0.4,RowSideColors = t(myCols), col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 0.9, KeyValueName = "Normalized Reads",Rowv = NA, Colv = NA, dendrogram = "none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(8,5))
  dev.off()
  
  ### perform clustering on both!
  
  mtrx = tmpMtrx
  mt = as.matrix(mtrx[,1:(ncol(mtrx)-1)])
  head(mt)
  
  ### generate colors
  myCols = getColors(mtrx, colNames = c("1","2","3","4","5","6","subplate","generalcortex"), coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black"))
  
  head(myCols)
  mtrx$LAYER
  nrow(mtrx)
  nrow(myCols)
  
  tmpMt = mt
  colnames(mt)
  
  pdf(paste(clstrFolder,"/Clustered_heatMap_allsamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4,RowSideColors = t(myCols), col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 0.9, KeyValueName = "Normalized Reads",Rowv = TRUE, Colv = TRUE, 
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(8,5))
  dev.off()
  
  colnames(mt)
  mt = mt[,-c(1,3,6)]
  colnames(mt)
  
  pdf(paste(clstrFolder,"/Clustered_heatMap_ControlSamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4,RowSideColors = t(myCols), col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 0.9, KeyValueName = "Normalized Reads",Rowv = TRUE, Colv = TRUE, 
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(8,5))
  dev.off()
  
  
  
  ##### Stem Cells
  
  if(!("mtrx.rda" %in% dir("rdaFiles/"))){
    mtrx = allNormData
    rownames(mtrx) = mtrx$Row.names
    mtrx$Row.names = NULL
    mtrx$ENSEMBL = NULL
    mtrx = mtrx[mtrx$SYMBOL %in% candidateGenes$SYMBOL,]
    mtrx = mtrx[!is.na(mtrx[,1]),]
    mtrx = merge(mtrx,candidateGenes,by = "SYMBOL", all = TRUE)
    mtrx = mtrx[!duplicated(mtrx$SYMBOL),]
    rownames(mtrx) = mtrx$SYMBOL
    mtrx$SYMBOL = NULL
    mtrx$ENSEMBL = NULL
    save(mtrx, file = "rdaFiles/mtrx.rda")
  }else{
    load("rdaFiles/mtrx.rda")
  }
  
  head(mtrx)
  
  mtrx = mtrx[mtrx$LAYER == "stemcell",]
  
  namesOrdered = c("SOX2","SOX9","GFAP","MSI1","PROM1","ASCL1")
  
  mtrx = mtrx[rownames(mtrx) %in% namesOrdered,]
  
  tmpMtrx = mtrx # keep for later use!
  
  tmpMat = mtrx
  tmpMat = tmpMat[rownames(tmpMat) %in% namesOrdered,]
  tmpMat = tmpMat[order(match(rownames(tmpMat),namesOrdered)),]
  
  mtrx  = tmpMat
  mt = as.matrix(mtrx[,1:(ncol(mtrx)-1)])
  head(mt)
  
  tmpMt = mt
  colnames(mt)
  
  pdf(paste(clstrFolder,"/stemcell_NoChangeInOrder_heatMap_allsamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4, col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "Normalized Reads",Rowv = NA, Colv = NA, dendrogram = "none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.0001,0.0001), lhei = c(0.0001,0.0001), margins = c(8,5))
  dev.off()
  
  colnames(mt)
  mt = mt[,-c(1,3,6)]
  colnames(mt)
  
  pdf(paste(clstrFolder,"/stemcell_NoChangeInOrder_heatMap_ControlSamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4, col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "Normalized Reads",Rowv = NA, Colv = NA, dendrogram = "none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.0001,0.0001), lhei = c(0.0001,0.0001), margins = c(8,5))
  dev.off()
  
  mtrx = tmpMtrx
  mt = as.matrix(mtrx[,1:(ncol(mtrx)-1)])
  head(mt)
  
  tmpMt = mt
  colnames(mt)
  
  pdf(paste(clstrFolder,"/stem_cell_Clustered_heatMap_allsamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4, col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 0.9, KeyValueName = "Normalized Reads",Rowv = TRUE, Colv = TRUE, 
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.0001,0.0001), lhei = c(0.0001,0.0001), margins = c(8,5))
  dev.off()
  
  colnames(mt)
  mt = mt[,-c(1,3,6)]
  colnames(mt)
  
  pdf(paste(clstrFolder,"/stemcell_Clustered_heatMap_ControlSamples_",rangeLim,".pdf", sep = ""))
  heatmap.3(mt, hclustfun = hclustfunc, distfun = distfunc, breaks = getBreaksForKey(maxLim = max(mt), rangeLim = rangeLim), RowSideColorsSize = 0.4, col = getBreaksColors(rangeLim),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 0.9, KeyValueName = "Normalized Reads",Rowv = TRUE, Colv = TRUE, 
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.0001,0.0001), lhei = c(0.0001,0.0001), margins = c(8,5))
  dev.off()
  
  #######
  
  if(!("allLogData.rda" %in% dir("rdaFiles/"))){
    log2Reads = rsDfr$dfr$deseq$outDfr[,1:(6*laCo)] 
    head(log2Reads)
    allLogData = merge(log2Reads,ann, by = "row.names", all = TRUE)
    save(allLogData, file = "rdaFiles/allLogData.rda")
  }else{
    load("rdaFiles/allLogData.rda")
  }
  
  head(allLogData)
  
  if(!("mtrxLog.rda" %in% dir("rdaFiles/"))){
    mtrxLog = allLogData
    rownames(mtrxLog) = mtrxLog$Row.names
    mtrxLog$Row.names = NULL
    mtrxLog$ENSEMBL = NULL
    mtrxLog = mtrxLog[mtrxLog$SYMBOL %in% candidateGenes$SYMBOL,]
    mtrxLog = mtrxLog[!is.na(mtrxLog[,1]),]
    head(mtrxLog)
    mtrxLog = merge(mtrxLog,candidateGenes,by = "SYMBOL", all = TRUE)
    mtrxLog = mtrxLog[!duplicated(mtrxLog$SYMBOL),]
    rownames(mtrxLog) = mtrxLog$SYMBOL
    mtrxLog$SYMBOL = NULL
    mtrxLog$ENSEMBL = NULL
    
    head(mtrxLog)
    
    save(mtrxLog, file = "rdaFiles/mtrxLog.rda")
  }else{
    load("rdaFiles/mtrxLog.rda")
  }
  
  head(mtrxLog)
  
  mtrxLog = mtrxLog[mtrxLog$LAYER != "stemcell",]
  
  namesOrdered = c("RELN","MARCKSL1","PLXND1","SATB2",
                   "MDGA1","POU3F2","TLE3","LHX2","POU3F3",
                   "TLE1","CUX1","MEF2C","DTX4","SLITRK1",
                   "UNC5D","CUX2","RORB","KCNIP2","CYP39A1","S100A10",
                   "FOXG1","LIX1L","OMA1","CNTN6","RAC3","FOXO1","OPN3","LDB2","CRYM","DKK3","TLE4","SEMA3E",
                   "OTX1","SOX5","LXN","FOXP2","CRIM1","PCP4","CTGF","DLX1","IGFBP4","ID2")
  
  mtrxLog = mtrxLog[rownames(mtrxLog) %in% namesOrdered,]
  
  tmpMtrxLog = mtrxLog # keep for later use!
  
  head(tmpMtrxLog)
  
  tmpMatLog = mtrxLog
  tmpMatLog = tmpMatLog[rownames(tmpMatLog) %in% namesOrdered,]
  tmpMatLog = tmpMatLog[order(match(rownames(tmpMatLog),namesOrdered)),]
  
  mtrxLog  = tmpMatLog
  colnames(mtrxLog)
  mt = as.matrix(mtrxLog[,seq(2, (ncol(mtrxLog)-1),6)])
  colnames(mt)
  colnames(mt) = gsub("\\.log2FoldChange$","",colnames(mt))
  head(mt)
  
  ### generate colors
  myCols = getColors(mtrxLog, colNames = c("1","2","3","4","5","6","subplate","generalcortex"), coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black"))
  
  head(myCols)
  mtrxLog$LAYER
  nrow(mtrxLog)
  nrow(myCols)
  
  tmpMt = mt
  head(mt)
  
  myBreaks <- c(seq(-3,-1,0.1),seq(-1,1,0.1),seq(1,3,0.1))
  myCol <- colorRampPalette(c("green", "yellow", "red"))(n = (length(myBreaks)-1))
  colnames(mt)
  
  pdf(paste(clstrFolder,"/log2Fold_NoChangeInOrder_heatMap_allsamples.pdf", sep = ""))
  heatmap.3(as.matrix(mt), hclustfun = hclustfunc, distfun = distfunc, breaks = myBreaks, col = myCol, RowSideColorsSize = 0.4, RowSideColors = t(myCols),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "log2FoldChange",Rowv = NA, Colv = NA, dendrogram = "none",symm=F,symkey=F,symbreaks=T, scale="none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(11,6))
  dev.off()
  
  pdf(paste(clstrFolder,"/log2Fold_Clustered_heatMap_allsamples.pdf", sep = ""))
  heatmap.3(as.matrix(mt), hclustfun = hclustfunc, distfun = distfunc, breaks = myBreaks, col = myCol, RowSideColorsSize = 0.4, RowSideColors = t(myCols),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "log2FoldChange",Rowv = TRUE, Colv = TRUE, symm=F,symkey=F,symbreaks=T, scale="none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(11,6))
  dev.off()
  
  
  colnames(mt)
  mt = mt[,c(3,4,5)]
  colnames(mt)
  
  pdf(paste(clstrFolder,"/log2Fold_NoChangeInOrder_ab42_VS_control_heatMap.pdf", sep = ""))
  heatmap.3(as.matrix(mt), hclustfun = hclustfunc, distfun = distfunc, breaks = myBreaks, col = myCol, RowSideColorsSize = 0.4, RowSideColors = t(myCols),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "log2FoldChange",Rowv = NA, Colv = NA, dendrogram = "none",symm=F,symkey=F,symbreaks=T, scale="none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(11,6))
  dev.off()
  
  pdf(paste(clstrFolder,"/log2Fold_Clustered_ab42_VS_control_heatMap.pdf", sep = ""))
  heatmap.3(as.matrix(mt), hclustfun = hclustfunc, distfun = distfunc, breaks = myBreaks, col = myCol, RowSideColorsSize = 0.4, RowSideColors = t(myCols),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "log2FoldChange",Rowv = TRUE, Colv = TRUE, symm=F,symkey=F,symbreaks=T, scale="none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(11,6))
  dev.off()
  
  if(!("mtrxLog.rda" %in% dir("rdaFiles/"))){
    mtrxLog = allLogData
    
    rownames(mtrxLog) = mtrxLog$Row.names
    mtrxLog$Row.names = NULL
    mtrxLog$ENSEMBL = NULL
    mtrxLog = mtrxLog[mtrxLog$SYMBOL %in% candidateGenes$SYMBOL,]
    mtrxLog = mtrxLog[!is.na(mtrxLog[,1]),]
    mtrxLog = merge(mtrxLog,candidateGenes,by = "SYMBOL", all = TRUE)
    mtrxLog = mtrxLog[!duplicated(mtrxLog$SYMBOL),]
    rownames(mtrxLog) = mtrxLog$SYMBOL
    mtrxLog$SYMBOL = NULL
    mtrxLog$ENSEMBL = NULL
    
    head(mtrxLog)
    
    save(mtrxLog, file = "rdaFiles/mtrxLog.rda")
  }else{
    load("rdaFiles/mtrxLog.rda")
  }
  
  head(mtrxLog)
  
  mtrxLog = mtrxLog[mtrxLog$LAYER == "stemcell",]
  
  namesOrdered = getNameOrders(mtrxLog)
  mtrxLog = mtrxLog[rownames(mtrxLog) %in% namesOrdered,]
  
  tmpMtrxLog = mtrxLog # keep for later use!
  
  head(tmpMtrxLog)
  
  tmpMatLog = mtrxLog
  tmpMatLog = tmpMatLog[rownames(tmpMatLog) %in% namesOrdered,]
  tmpMatLog = tmpMatLog[order(match(rownames(tmpMatLog),namesOrdered)),]
  
  mtrxLog  = tmpMatLog
  
  mt = as.matrix(mtrxLog[,seq(2, (ncol(mtrxLog)-1),6)])
  mt
  colnames(mt) = gsub("\\.log2FoldChange$","",colnames(mt))
  head(mt)
  ### generate colors
  myCols = getColors(mtrxLog, colNames = c("1","2","3","4","5","6","subplate","generalcortex"), coloRs = c("skyblue","brown","red","orange","pink","violet","purple","gray","green","blue","yellow","black"))
  
  head(myCols)
  mtrxLog$LAYER
  nrow(mtrxLog)
  nrow(myCols)
  
  tmpMt = mt
  max(mt, na.rm = TRUE)
  head(mt)
  
  colnames(mt)
  
  pdf(paste(clstrFolder,"/log2Fold_stemcell_NoChangInOrder_heatMap_allsamples.pdf", sep = ""))
  heatmap.3(as.matrix(mt), hclustfun = hclustfunc, distfun = distfunc, breaks = myBreaks, col = myCol, RowSideColorsSize = 0.4, RowSideColors = t(myCols),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "log2FoldChange",Rowv = NA, Colv = NA, dendrogram = "none",symm=F,symkey=F,symbreaks=T, scale="none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(11,6))
  dev.off()
  
  colnames(mt)
  
  pdf(paste(clstrFolder,"/log2Fold_stemcell_Clustered_heatMap_allsamples.pdf", sep = ""))
  heatmap.3(as.matrix(mt), hclustfun = hclustfunc, distfun = distfunc, breaks = myBreaks, col = myCol, RowSideColorsSize = 0.4, RowSideColors = t(myCols),
            key = F,rowsep = 1:(nrow(mt)+1),trace = "none",  colsep=1:(ncol(mt)+1),sepcolor='white',density.info='none',sepwidth = c(0.0001,0.0001),
            keysize = 1, KeyValueName = "log2FoldChange",Rowv = TRUE, Colv = TRUE, dendrogram = "none",symm=F,symkey=F,symbreaks=T, scale="none",
            cexRow = 0.7, cexCol = 0.7, lwid = c(0.2,0.2), margins = c(11,6))
  dev.off()
}  


