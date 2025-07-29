#########################################################################
#### Script to download and clean GWAS for heterogeneity comparisons ####


#### Steps
      
    # 1. Read in data
          # read in source list
          # Access IEU Open GWAS
              # Download associations
                # Do for GWAS 1
                # Get clumped Top hits from both
                # Extract subset of SNPs in GWAS 1 from GWAS 2
          # Download EBI GWAS

    # 2. data cleaning
          


######################################################################



######## Packages ######## 

  library(data.table)
  library(dplyr)
  library(metafor)
  library(ieugwasr)
  library(gwasglue)
  library(httr)
  library(curl)
  library(gwasrapidd)
  library(VariantAnnotation)


##########################

setwd("../../data/")
token <- readLines("token")[1]
timeout(560) # set timeout for downlaods
###########################################################

  ################## 1. Read in data ###################

###########################################################

latAm_gwas <- fread("multi_gwas_across_sources.txt")

euro_gwas <- fread("euro_gwas_pairs.txt")

IEU_ids <- c(latAm_gwas[latAm_gwas$source == "IEU",]$id, euro_gwas$id)

EBI_ids <- latAm_gwas[latAm_gwas$source == "EBI",]$id


################## load IEU GWAS
no_vcf_list <- list()
for (id in IEU_ids){
  
  if (! file.exists(paste0("vcfs/", id, ".vcf.gz"))) {
    
    url <- paste0("https://gwas.mrcieu.ac.uk/files/",id,"/",id,".vcf.gz")
    
    ## check URl exists
    res <- HEAD(url)
    
    if (status_code(res) == 200) {
      message(paste(id, "URL exists."))
      curl_download(url, destfile = paste0("vcfs/",id, ".vcf.gz"))
    } else {
      message(paste(id, "URL does not exist or is not reachable."))
      no_vcf_list <- append(no_vcf_list, id)
    }
    
    
    
  }
  
  
}

### add missing vcf files to EBI ids to check there


########### Need to extract all studies for range required in GWAS catalog 

ebi_pmids  <- as.character(latAm_gwas[latAm_gwas$source == "EBI" | latAm_gwas$id %in% no_vcf_list,]$pmid)

ebi_full_study_list <- get_studies(pubmed_id  = ebi_pmids)


EBI_ids <- c(EBI_ids, no_vcf_list)


EBI_studies <- latAm_gwas[latAm_gwas$id %in% EBI_ids,]
EBI_studies$id <- gsub("ebi-a-", "", EBI_studies$id)

EBI_data_locs <- fread("gwas_cat_locs.txt", header=F)

EBI_studies$data_path <- sapply(EBI_studies$id, function(study_id) {
  match_idx <- grep(study_id, EBI_data_locs$V1, fixed = TRUE)
  if (length(match_idx) > 0) {
    return(EBI_data_locs$V1[match_idx[1]])  # Return the first match (should be unique)
  } else {
    return(NA)  # No match foundS
  }
})


EBI_studies <-EBI_studies[EBI_studies$trait != "asthma",]  # remove asthma as ancestry was miss-labelled



#### loop download EBI GWAS
for (study_id in EBI_studies$id){
  
  if (! file.exists(paste0("vcfs/", study_id, ".vcf.gz"))) {
    
    study <- EBI_studies[EBI_studies$id == study_id , ]
    url <- paste0("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", study$data_path)

    ## check URl exists
    res <- HEAD(url)
    
    if (status_code(res) == 200) {
      message(paste(study_id, "URL exists."))
      curl_download(url, destfile = paste0("vcfs/",study_id, ".vcf.gz"))
    } else {
        message(paste(study_id, "URL does not exist"))
    }
    
    
    
  }
  
  
}

###########################################################

 ################## 1. data cleaning ###################

###########################################################

t <- readVcf("vcfs/ebi-a-GCST004605.vcf.gz")

files <- list.files("vcfs/", full.names = TRUE)




###########################################################

################ 1. heterogeneity testing #################

###########################################################


for (current_trait in unique_traits)  {
  
  g_l <- latAm_gwas[latAm_gwas$trait == current_trait,]
  g_l <- g_l[order(g_l$sample_size),]
  
  GWAS_1 <- readVcf(paste0("vcfs/",g_l[1,]$id, ".vcf.gz"))
  GWAS_2 <- readVcf(paste0("vcfs/",g_l[2,]$id, ".vcf.gz"))
  
  if (g_l[1,]$source == "IEU"){
    GWAS_1 <- GWAS_1[COL_SELECTIOON]
  
  }
  
}

