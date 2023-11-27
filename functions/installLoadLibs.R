#updated 26.08.2017
installLoadLibs <- function(){
	ReqLibs = c("AnnotationForge","BiocParallel","biomaRt","clusterProfiler","calibrate","data.table",
	"DESeq2","devtools", "gplots","plyr","dtplyr","fdrtool","genefilter", #"FactoMineR",
	"geneLenDataBase","geneplotter","GenomicFeatures","ggplot2","ggrepel","GO.db",
	"goseq","GOstats","gplots","gdata","GSEABase","KEGG.db", "KEGGREST",     
	"LSD", "parallel","pheatmap","plyr","pathview","PerformanceAnalytics",
	"plotrix","png","RColorBrewer","ReportingTools","topGO","reshape2",
	"Rgraphviz","splitstackshape","scales","statmod","stringr","UpSetR","xlsx","openxlsx","XML",	
	"xtable")
	if(!is.element("BiocInstaller",installed.packages()[,1])){
		source("http://www.bioconductor.org/biocLite.R")
		biocLite()
		library("BiocInstaller")
	}else{
		library("BiocInstaller")
	}
	for (i in 1:length(ReqLibs)){
		ReqLib = ReqLibs[i]
		if(!is.element(ReqLib, installed.packages()[,1])){
			biocLite(ReqLib)
		}
	}
	lapply(ReqLibs, require, character.only = T)
}
