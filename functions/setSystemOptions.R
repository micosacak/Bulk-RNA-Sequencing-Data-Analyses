#updated on 20.10.2017
setSystemOptions <- function(min8 = 2){
	if(Sys.info()["sysname"][[1]] == "Linux") {
		endsWith <- stringi::stri_endswith # endsWith does not work on Linux!!! instead, use stri_endswith
		register(MulticoreParam(detectCores()-min8))
	}else if(Sys.info()["sysname"][[1]] == "Windows"){
		register(SnowParam(detectCores()-min8))
	}
}
