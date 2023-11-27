#updated on 26.08.2017
getOrgInfo <- function(orgId){
	library(paste0("org.",orgId,".eg.db"),character.only = TRUE)
	results <- list("orgId" = NULL, "org3L" = NULL, "orgDb" = NULL, "orgNm" = NULL, "orgSp" = NULL, "orgGb" = NULL, "orgEs" = NULL)
	if(orgId == "Dr"){
		results$orgEs <- "drerio_gene_ensembl"  # ensembl dataset name
		results$orgId <- "Dr"                   # 2-letters organism ID
		results$org3L <- "dre"                  # 3-letters organism ID
		results$orgDb <- "org.Dr.eg.db"         # orgnism annotaion DB, as org.Xx.eg.db
		results$orgNm <- "drerio"               # organism name shortened
		results$orgSp <- "Danio rerio"          # organism names as species full names!
		results$orgGb <- "danRer10"             # organism genome build
		library("org.Dr.eg.db")
	}else if(orgId == "Hs"){
		results$orgEs <- "hsapiens_gene_ensembl"
		results$orgId <- "Hs"
		results$org3L <- "hsa"
		results$orgDb <- "org.Hs.eg.db"
		results$orgNm <- "hsapiens"
		results$orgSp <- "Homo sapeins"
		results$orgGb <- "hg19"
		library("org.Hs.eg.db")
	}else if(orgId == "Mm"){
		results$orgEs <- "musculus_gene_ensembl"
		results$orgId <- "Mm"
		results$org3L <- "mmu"
		results$orgDb <- "org.Mm.eg.db"
		results$orgNm <- "mmusculus"
		results$orgSp <- "Mus musculus"
		results$orgGb <- "mmu10"
		library("org.Mm.eg.db")
	}else if(orgId == "Rn"){
		results$orgEs <- "rnorvegicus_gene_ensembl"
		results$orgId <- "Rn"
		results$org3L <- "rnu"
		results$orgDb <- "org.Rn.eg.db"
		results$orgNm <- "rnorvegicus"
		results$orgSp <- "Rattus norvegicus"
		results$orgGb <- "rn6"
		library("org.Rn.eg.db")
	}else{
		stop(paste0("The orgId = ",orgId," does not have 3-letters, annotaion databass, ...\n",
		"please put these names in the function getOrgInfo using Dr, Hs, Mm and Rn as examples!!!"
		))
	}
	return(results)
}
