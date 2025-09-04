############################################################
#### Script to perform data manipulation  for BMI data  ####


#### Steps
# read in trait list
# split into LA and Eur studies
# read in BMI GWAS and format iteratively
# Harmonise phenotypes
# get clumped top hits fromall GWAS, make master SNP list


###########################################################

library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(genetics.binaRies, lib.loc= "/user/home/kb22541/R/x86_64-conda-linux-gnu-library/4.4")
library(gwasvcf)
library(VariantAnnotation)


setwd("/user/work/kb22541/local_anc/data/")
token <- readLines("token")[1]



########################

  # Define column names

########################

gwas_cat_cols <- c("SNP","POS","CHROM","REF","ALT","beta","se","p_value","RAF")
ieu_cols <- c("rsid","start","seqnames","REF","ALT","ES","SE","LP","AF")
output_cols <- c("MARKERNAME","POSITION","CHROMOSOME","EA","NEA","BETA","SE","PVAL","EAF","N")



trait_list <- fread("multi_gwas_across_sources.txt")
trait_list <- trait_list[(trait_list$trait %in%  c("body mass index", "Body mass index")) & trait_list$sample_size > 5000, ]

trait_list_LA <- trait_list[trait_list$ancestry == "Hispanic or Latin American",]
trait_list_EUR <- trait_list[trait_list$ancestry == "European",]

## Check if files exist


# results <-c()
# 
# ## get top hits list
# for (file in trait_list$id){
#   
# 
#   filepath <- (paste0("vcfs/", file, ".vcf.gz"))
#   
#   if (file.exists(filepath)){
#     
#     message( "File ", file, " exists: TRUE")
#     
#     if (trait_list[trait_list$id == file,]$source == "IEU"){
#     
#       tophits <-suppressMessages(vcf_to_tibble(query_gwas(filepath, pval=5e-7)))
#       names(tophits)[names(tophits) == "rsid"] <- "SNP"
#       
#       
#     } else {
#       
#       tophits <- fread(filepath)   # reads tab-separated, gzipped
#       tophits <- tophits[tophits$p_value < 5e-7, ]  
#       
#     }
#     
#     
#     
#     if (length(tophits$SNP > 1)){
#       
#       message(length(tophits$SNP), "top hits for file ", file)
#       
#       results <-c(results, tophits$SNP)
#     } else{message("less than 1 top hit for ", file)}
#     
#     
#   } else{message(paste("File", file, "does not exist.")) }
#   
# }
# 
# results <- unique(results)
# message("Total SNPs = ", length(results) ," \nFinsihed tophits extraction, making subsets..")


for (file in trait_list$id){
  
     filepath <- (paste0("vcfs/", file, ".vcf.gz"))
  
     
     if (trait_list[trait_list$id == file,]$source == "EBI"){
        cols <- gwas_cat_cols
        main <- fread(filepath)
        print(paste(file, "columns: ", colnames(main)))

        main <- main[, cols, with = FALSE]

        subset <- main[main$p_value <=5e-7,]
     }
        
        
      
    
     else { 
        cols <- ieu_cols

        
        main <-(vcf_to_tibble(readVcf(filepath)))
        print(paste(file, "columns: ", colnames(main)))
        
        main <- main[, cols, with = FALSE]  
        # subset <- vcf_to_tibble(query_gwas(main, pval = 5e-7))
     }
    
     
    # message(length(subset$rsid), "of",length(results), "extracted for ", file)
    # 
    # subset <- subset[,..cols]
    # subset$N <- 1000
    
    main$N <- max(trait_list[trait_list$id == file,]$sample_size)
    
    colnames(main) <- output_cols
    write.table(main, paste0("subset/", file,"_main.txt"), sep= "\t", row.names=F, quote=F)
    # colnames(subset) <- output_cols
    
    # write.table(subset, paste0("subset/", file, ".txt"), sep = "\t", row.names=F, quote=F)
    
}







