###########################################################

#### Script to perform GWAS pair heterogeneity testing ####


#### Steps
      
    # 1. Read in data
          # read in source list
          # Access IEU Open GWAS
              # Download associations
                # Do for GWAS 1
                # Get clumped Top hits from both
                # Extract subset of SNPs in GWAS 1 from GWAS 2
          # Download EBI GWAS

        
    # 2. calculate heterogeneity stats
          # Chochrane's Q
          # I squared
          # Replication rates

###########################################################



######## Packages ######## 

  library(data.table)
  library(dplyr)
  library(metafor)
  library(ieugwasr)
  library(gwasglue)

##########################



###########################################################

  ################## 1. Read in data ###################

###########################################################

latAm_gwas <- fread("multi_gwas_across_sources.txt")

euro_gwas <- fread("euro_gwas_pairs.txt")

IEU_ids <- c(latAm_gwas[latAm_gwas$source == "IEU",]$id, euro_gwas$id)

EBI_ids <- latAm_gwas[latAm_gwas$source == "EBI",]$id


################## load ieu GWAS














