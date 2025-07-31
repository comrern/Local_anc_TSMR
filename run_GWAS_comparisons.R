#########################################################################
#### Script to perform heterogeneity calculations and fine mapping  ####


#### Steps
  # read in trait list
    # Iterate through all traits
      # Harmonise phenotypes
      # get clumped top hits from GWAS 1, het calcs against same SNPs in GWAS 2
      # Do the above in reverse, compare top hits in GWAS 2 to 1
      # Calculate observed / expected replication rate
      # Repeat the above for fine Mapped SNPs (CAMERa)


######################################################################

library(CAMERA)
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(gwasvcf)
library(VariantAnnotation)

setwd("../../data/")
token <- readLines("token")[1]


###########################################################

    ################## Functions ##################

###########################################################

pheno_harm <- function(ids){
  
    o <- lapply(ids, function(i) {
      tryCatch(
        {
          suppressMessages(exp <- unique(TwoSampleMR::extract_instruments(outcomes = i ,opengwas_jwt = token)))
          if (length(exp$SNP) < 2){
            message("Not enough significant SNPs for MR, reducing threshold")
            suppressMessages(exp <- unique(TwoSampleMR::extract_instruments(outcomes = i , p1 = 5e-7, opengwas_jwt = token)))
          }
          
          other_ids <- ids[!ids %in% i]
          
          o <- lapply(other_ids, function(j) {
            suppressMessages(out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP, outcomes = j, opengwas_jwt = token))
            suppressMessages(d <- TwoSampleMR::harmonise_data(exp, out))
            
            res <- suppressMessages(TwoSampleMR::mr(d, method = "mr_ivw")) %>%
              {
                dplyr::tibble(Reference = i, Replication = j, nsnp = .$nsnp, agreement = .$b, se = .$se, pval = .$pval)
              }
            
            message(paste0("Instrument associations between ", i, " and ", j, " is ", round(res$agreement, 3), "; NSNP=", res$nsnp))

            return(res)
          })
          return(o %>% dplyr::bind_rows())
        },
        error = function(e) {
          cat("Unable to perform phenotype harmonisation :", conditionMessage(e), "\n")
        }
      )
    })
    o <- (o %>% dplyr::bind_rows())

    return(o)  
  
}

heterogeneity <- function(data){
  
  Q_df <-  data.frame(ID = numeric(),
                      Qsnp = numeric(),
                      Qp = numeric(),
                      X1_b = numeric(),
                      X2_b = numeric(),
                      stringsAsFactors = FALSE)
  
  for (i in 1:nrow(data)){
    
    betas <- data[i,c("beta.x","beta.y")] 
    ses <- data[i,c("beta.x","beta.y")]
    
    w <- 1 / (ses)^2 # get weights
    ivw_b <- sum(betas * w) / sum(w) # ivw betas
    se <- sqrt(1 / sum(w)) # ivw se
    
    Q <- sum(w * (betas - ivw_b)^2)
    
    
    df <- length(betas) -1
    
    Qpval <- stats::pchisq(Q, df, lower.tail=FALSE)
    
    Q_df[i, ] <- list(i, Q, Qpval, betas[1, 1], betas[1, 2])
    
    
    
  }
  
  Qdf <- (length(Q_df$ID)  - 1)
  Q_all <- sum(Q_df$Qsnp)
  
  
  Q_total <-  c("Qsum",
                sum(Q_df$Qsnp), 
                pchisq(Q_all, df =  Qdf, lower.tail = FALSE),
                (mean(Q_df$Qp < 0.05) * 100),
                max(0, ((Q_all - Qdf)/ Q_all) * 100)
  )
  
  return(Q_total)  
  
  
}


get_instruments <- function(ids_f){
  
  g1 <- readVcf(paste0("vcfs/",ids_f[1], ".vcf.gz"))
  g2 <- readVcf(paste0("vcfs/",ids_f[2], ".vcf.gz"))
  
  g1_tophits <- vcf_to_tibble(query_gwas(g1, pval = 5e-8))
  g2_tophits <- vcf_to_tibble(query_gwas(g2, pval = 5e-8))
  
  g2_g1tophits <- vcf_to_tibble(query_gwas(g2, rsid = g1_tophits$rsid))
  g1_g2tophits <- vcf_to_tibble(query_gwas(g1, rsid = g2_tophits$rsid))
  
  g1_merge <- rbind(g1_tophits, g1_g2tophits[g1_g2tophits$rsid %in% g1_tophits$rsid,])
  g2_merge <- rbind(g2_tophits, g2_g1tophits[g2_g1tophits$rsid %in% g2_g1tophits$rsid,])
  
  
  ##### extract regions for tophit SNPs in either gwas
  regions <- paste0(g1_merge$seqnames, ":", g1_merge$start - 50000, "-", g1_merge$end + 50000)
  
  regions <- lapply(regions, function(r){
    
      a <- rbind(vcf_to_tibble(query_gwas(g1, chrompos = r)),
                 vcf_to_tibble(query_gwas(g2, chrompos = r)))  %>%
          dplyr::arrange(start) %>%
          dplyr::bind_rows()
        message(nrow(a))
   
      a <- lapply(ids_f, function(i) {
        subset(a, id == i) %>%
          dplyr::filter(!duplicated(rsid))
      })
      rsids <- Reduce(intersect, lapply(a, function(x) x$rsid))
      a <- lapply(a, function(x) {
        subset(x, rsid %in% rsids)
      })
      ALT <- a[[1]]$ALT
      a <- lapply(a, function(x) {
        index <- x$ALT != ALT
        if (sum(index) > 0) {
          x$ES[index] <- x$ES[index] * -1
          REF <- x$REF[index]
          x$REF[index] <- x$ALT[index]
          x$ALT[index] <- x$REF[index]
          x$AF[index] <- 1 - x$AF[index]
          x <- subset(x, REF == a[[1]]$REF)
        }
        return(x)
      })
      rsids <- Reduce(intersect, lapply(a, function(x) x$rsid))
      a <- lapply(a, function(x) {
        subset(x, rsid %in% rsids) %>%
          dplyr::arrange(seqnames, start)
      })
      return(a)
      })
  
  ## merge all regions DFs
  
  # Extract all the first data frames (i.e., [[1]] from each sublist)
  r1 <- lapply(regions, function(x) x[[1]])
  
  # Extract all the second data frames (i.e., [[2]] from each sublist)
  r2 <- lapply(regions, function(x) x[[2]])
  
  # Combine them
  r1_combined <- do.call(rbind, r1)
  r2_combined <- do.call(rbind, r2)
  
  return(list(g1_raw = g1_merge, g2_raw = g2_merge, r1_fema = r1_combined, r2_fema = r2_combined ))
  
}

run_fema <- function(betas, ses) {
  w <- 1 / ses^2
  beta <- rowSums(betas * w) / rowSums(w, na.rm=TRUE)
  se <- sqrt(1 / rowSums(w, na.rm=TRUE))
  z <- abs(beta / se)
  p <- pnorm(z, lower.tail = FALSE)
  nstudy <- apply(betas, 1, \(x) sum(!is.na(x)))
  return(tibble(nstudy, p, z=z))
}
  
  
  
}

heterogeneity_calcs <- function(df1, df2, method){
  
  if (method == "raw"){
      
      g1xg2 <- merge(df1[df1$LP > 8,], df2, by= "rsid")
      g2xg1 <- merge(df2[df2$LP > 8,], df1, by= "rsid")

      g1xg2 <- heterogeneity(g1_merge)
      g2xg1 <- heterogeneity(g2_merge)
        
      out <- data.frame(rbind(g1xg2, g2xg1))
      out[,1] <- ids
      out$method = "raw"
      
  }    
  
  
  if (method == "fema"){
    
    fema_dat_merge <- merge(df1, df2, by = "rsid")
    
    FEMA_SNPS <- run_fema( as.matrix(fema_dat_merge[,"beta.x","beta.y"], fema_dat_merge["se.x","se.y"]) )
    
    out <- heterogeneity(FEMA_SNPS)
    out$method = "fema"                       
  }
  
  return(df)
  
}



###########################################################

################## 1. Read in trait list ##################

###########################################################


latAm_gwas <- fread("multi_gwas_across_sources.txt")
euro_gwas <- unique(fread("euro_gwas_pairs.txt"))
euro_gwas <- euro_gwas[, c("id","pmid","author","trait","population","sample_size")]
euro_gwas$source <- "IEU"
colnames(euro_gwas)[5] <- "ancestry"

# latAm_gwas <- latAm_gwas[latAm_gwas$id != "ebi-a-GCST008053",]


###########################################################

################ 1. heterogeneity testing #################

###########################################################
## function structure
# Need to separate IEU data (can use remote) from open GWAS

latAm_gwas <- latAm_gwas[latAm_gwas$source == "IEU" & latAm_gwas$trait != "asthma",] # remove asthma as this was matched to GWAS catalogue 
unique_traits <- unique(latAm_gwas$trait)

all_phen <- data.frame()
het <- data.frame()
for (current_trait in unique_traits)  {
  
  g_l <- latAm_gwas[latAm_gwas$trait == current_trait,]
  g_l <- g_l[order(g_l$sample_size, decreasing = T),]
  
  ids_m <-c(g_l[1,]$id, g_l[2,]$id)
  
  
  # p <- pheno_harm(ids_m)
  # all_phen <- dplyr::bind_rows(all_phen, p)
  
  instruments <- get_instruments(ids_m)
  
  
  h <- heterogeneity_calcs(instruments[g1_raw], instruments[g2_raw], "raw")
  het <- dplyr::bind_rows(het, h)
  
  h <- heterogeneity_calcs(instruments[r1_fema], instruments[r2_fema], "fema")
  het <- dplyr::bind_rows(het, h)
  
  
}









