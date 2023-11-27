drawPie <- function(Df = NULL, outputName = NULL, output_folder = NULL, pvalCol = NULL, labCol = NULL, countCol = NULL, totalCol = NULL, maxRes = 50, pvalCutoff = 0.05, ...){
  # draws a pie for gene ontology!
  # all fields are required!
  # Df         : data frame of the any result type!
  # outputName : the output folder name.
  # pvalCutoff : the pval cut off (default: 0.05)
  # maxRes     : maximum results to show on pie chart (default: nrow(Df)).
  # pvalCol    : which column has the p-values (or frd, etc.)
  # labCol     : which column has the labels.
  # countCol   : which column has the count values to draw pie chart.
  Df <- Df[Df[,pvalCol] < pvalCutoff,]
  Df <- Df[Df[,countCol] != 0,]
  Df <- Df[order(Df[,countCol], decreasing = TRUE),]
  Df[,labCol] <- strtrim(Df[,labCol],50)
  if(nrow(Df) > 0){
    tmpDf <- Df[order(Df[,pvalCol]),]
    if(nrow(tmpDf) > maxRes) tmpDf = tmpDf[1:maxRes,]
    pdf(paste(output_folder,paste("pVal_",outputName, sep = ""), sep = "/"))
    par(mar = c(5, 15, 4, 5) + 0.1)
    barplot(-log10(tmpDf[,pvalCol]), xlab = "-log10(pvalue)",horiz = TRUE, names.arg = tmpDf[,labCol], las = 1, cex.axis = 1, cex.names = 0.5)
    dev.off()
    
    if(nrow(Df) > maxRes){
      Df <- Df[1:maxRes,]
    }
    
    pdf(paste(output_folder,outputName, sep = "/"), width = 35, height = 28, pointsize = 14)
    pie(Df[,countCol], labels = paste(Df[,labCol],paste(Df[,countCol],Df[,totalCol], sep = "/"), sep = ":"), main = gsub(".pdf$","",outputName), radius = 0.5)
    dev.off()
  }
}