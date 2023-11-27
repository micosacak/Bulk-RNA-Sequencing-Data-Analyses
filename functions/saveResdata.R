#updated on 20.10.2017
saveResdata <- function(dds_res,dds,fdr_threshold, output_folder, output_file_name){
	resdata1 <- merge(as.data.frame(dds_res), as.data.frame(counts(dds, normalized = FALSE)), by = "row.names", sort = FALSE)
	names(resdata1)[1] <- "GeneID"
	resdata1_df <- resdata1[,2:ncol(resdata1)]
	rownames(resdata1_df) <- resdata1$GeneID
	resdata2 <- merge(resdata1_df, as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
	names(resdata2)[1] <- "GeneID"
	resdata <- resdata2[,2:ncol(resdata2)]
	rownames(resdata) <- resdata2$GeneID
	return(resdata)
}