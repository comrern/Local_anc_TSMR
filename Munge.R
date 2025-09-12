library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(data.table)
setwd("C:/Users/kb22541/Desktop/Analyses/Local_anc/data/mr_mega_dat/")




files <- list.files(pattern = "*\\.tsv$", full.names = TRUE)


for (file in files) {
  
    data <- fread(file)
    
    
    pval_cols <- c("Score.pval", "p", "PVAL", "pval")
    
    # find which are actually present
    present <- names(data) %in% pval_cols
    
    if (any(present)) {
      # rename all that match to "P"
      names(data)[present] <- "P"
    } else {
      stop("No recognised p-value column found")
    }
    
    data <- data[,c("MARKERNAME","CHROMOSOME","POSITION","NEA","EA","BETA","EAF","P","SE","N")]
    
    # colnames(data)["Score.pval"] <- "P"
    
    MungeSumstats::format_sumstats(data,ref_genome = "GRCh38" ,save_path = paste0(file, "_munged.txt"))


}


