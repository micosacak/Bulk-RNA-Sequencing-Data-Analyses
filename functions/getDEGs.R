getDEGs = function(dfrs, rs = "all", output_folder = "",lfc_threshold = 1.0, fdr_threshold = 0.1){
  res = list("up" = c(),"down" = c(),"all" = c())
  res["all"][[1]] = rownames(subset(dfrs, abs(log2FoldChange) >= lfc_threshold & padj < fdr_threshold))
  res["up"][[1]] = rownames(subset(dfrs, (log2FoldChange) <= -lfc_threshold & padj < fdr_threshold))
  res["down"][[1]] = rownames(subset(dfrs, (log2FoldChange) >= lfc_threshold & padj < fdr_threshold))
  return(res[rs][[1]])
}
