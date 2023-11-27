#updated on 26.08.2017
mergeAllDfr <- function(theListofDfrs){
	rs <- theListofDfrs[[1]]
	if(length(theListofDfrs) >= 2){
		for(i in 2:length(theListofDfrs)){
			rs <- merge(rs, theListofDfrs[[i]], by = "row.names", all.x = TRUE, all.y = TRUE)
			rownames(rs) <- rs$Row.names
			rs$Row.names <- NULL
		}
	}
	return(rs)
}