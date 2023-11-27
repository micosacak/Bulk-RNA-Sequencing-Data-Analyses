#updated on 20.10.2017
readInputFiles <- function(input_file = NULL, rmSuf = ".fastq.gz.sam", rmCol = -c(1:5), rowID = "Geneid"){
  rawCounts = read.delim(input_file, skip = 1, row.names = rowID)
  rawCounts <- rawCounts[,rmCol]
  colnames(rawCounts) <- gsub(rmSuf,"",colnames(rawCounts))
  colnames(rawCounts)
  return(rawCounts)
}