### draw venn schema

drawVenn <- function(dds = NULL, aCo = NULL, output_folder = "",lfc_threshold = 1.0, fdr_threshold = 0.1){
  ddsResults = vector("list", length = 3)
  names(ddsResults) = c("all","down","up")
  for(adu in c("all","down","up")){
    ddsResults[[adu]] = vector("list", length = length(aCo))
    names(ddsResults[[adu]]) = aCo
  }
  for(comp in aCo){
    trt_ctrl = unlist(strsplit(comp, split = "_vs_"))
    ddsRes = results(dds, contrast = c("condition", trt_ctrl))
    for(adu in c("all","up","down")) ddsResults[[adu]][[comp]] = getDEGs(as.data.frame(ddsRes), rs = adu)
  }  
  
  delistALL = ddsResults$all
  delistDW = ddsResults$down
  delistUP = ddsResults$up
  
  if(length(aCo) > 3){
    pdf(paste(output_folder,"venn_all_comparisons_ALL_degs.pdf", sep = "/"))
    DE_gns = UpSetR::fromList(delistALL)
    UpSetR::upset(DE_gns, order.by = "freq")
    dev.off()
    
    pdf(paste(output_folder,"venn_all_comparisons_UP_degs.pdf", sep = "/"))
    DE_gns = UpSetR::fromList(delistUP)
    UpSetR::upset(DE_gns, order.by = "freq")
    dev.off()
    
    pdf(paste(output_folder,"venn_all_comparisons_DOWN_degs.pdf", sep = "/"))
    DE_gns = UpSetR::fromList(delistDW)
    UpSetR::upset(DE_gns, order.by = "freq")
    dev.off()
  }else{
    pdf(paste(output_folder,"venn_all_comparisons_ALL_degs.pdf", sep = "/"))
    DE_gns = UpSetR::fromList(delistALL)
    UpSetR::upset(DE_gns, order.by = "freq")
    dev.off()
    
    pdf(paste(output_folder,"venn_3_comparisons_ALL_degs.pdf", sep = "/"))
    venn(delistALL)
    dev.off()
  
    
    pdf(paste(output_folder,"venn_all_comparisons_UP_degs.pdf", sep = "/"))
    DE_gns = UpSetR::fromList(delistUP)
    UpSetR::upset(DE_gns, order.by = "freq")
    dev.off()

    pdf(paste(output_folder,"venn_3_comparisons_DOWN_degs.pdf", sep = "/"))
    venn(delistDW)
    dev.off()
        
    pdf(paste(output_folder,"venn_all_comparisons_DOWN_degs.pdf", sep = "/"))
    DE_gns = UpSetR::fromList(delistDW)
    UpSetR::upset(DE_gns, order.by = "freq")
    dev.off()
    
    pdf(paste(output_folder,"venn_3_comparisons_UP_degs.pdf", sep = "/"))
    venn(delistUP)
    dev.off()
    
  }
  
}


