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

library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(gwasvcf)
library(CAMeRa, lib.loc = "/user/home/kb22541/R/x86_64-conda-linux-gnu-library/4.4/")
library(genetics.binaRies, lib.loc= "/user/home/kb22541/R/x86_64-conda-linux-gnu-library/4.4")
library(VariantAnnotation)

setwd("/user/work/kb22541/local_anc/data/")
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
    
    betas <- data[i,c("ES.x","ES.y")] 
    ses <- data[i,c("SE.x","SE.y")]
    
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

  
  g1_path <- paste0("vcfs/",ids_f[1], ".vcf.gz")
  g2_path <- paste0("vcfs/",ids_f[2], ".vcf.gz")
  
  # get top this
  tryCatch({
    
    print(paste("Reading", g1_path))
    g1_tophits <- suppressMessages(vcf_to_tibble(query_gwas(g1_path, pval=5e-8)))
    print(paste("Done, Reading", g2_path))
    g2_tophits <- suppressMessages(vcf_to_tibble(query_gwas(g2_path, pval=5e-8)))
    
    if (length(g1_tophits$rsid )< 1){
      g1_tophits <- vcf_to_tibble(query_gwas(g1_path, pval=5e-7))
    }
    if (length(g2_tophits$rsid) < 1){
      g2_tophits <- vcf_to_tibble(query_gwas(g2_path, pval=5e-7))
    }

    g1_tophits$pval <- 10^(-g1_tophits$LP)
    g2_tophits$pval <- 10^(-g2_tophits$LP)
    
    
    ### clump hits
    clump_list <- ld_clump(
      tibble(rsid = g1_tophits$rsid , pval = g1_tophits$pval ),
      plink_bin = get_plink_binary(),
      bfile = "AMR"
    )
    g1_tophits <- g1_tophits[g1_tophits$rsid %in% clump_list$rsid, ]
    
    clump_list <- ld_clump(
      tibble(rsid = g2_tophits$rsid , pval = g2_tophits$pval ),
      plink_bin = get_plink_binary(),
      bfile = "AMR"
    )
    g2_tophits <- g2_tophits[g2_tophits$rsid %in% clump_list$rsid, ]
    
    
    # Step 2: Keep only rows where chr is between 1 and 22
    g1_tophits <- g1_tophits[as.numeric(as.character(g1_tophits$seqnames)) %in% 1:22, ]
    g2_tophits <- g2_tophits[as.numeric(as.character(g2_tophits$seqnames)) %in% 1:22, ]
    
    
    
    g1_tophits$chrpos <- paste0(g1_tophits$seqnames, ":", g1_tophits$start, "-", g1_tophits$end )
    g2_tophits$chrpos <- paste0(g2_tophits$seqnames, ":", g2_tophits$start, "-", g2_tophits$end )
    
    g2_g1tophits <- vcf_to_tibble(query_gwas(g2_path, chrompos = g1_tophits$chrpos))
    g1_g2tophits <- vcf_to_tibble(query_gwas(g1_path, chrompos = g2_tophits$chrpos))
    
    g2_g1tophits$chrpos <- paste0(g2_g1tophits$seqnames, ":", g2_g1tophits$start, "-", g2_g1tophits$end )
    g1_g2tophits$chrpos <- paste0(g1_g2tophits$seqnames, ":", g1_g2tophits$start, "-", g1_g2tophits$end )
    
    
    g1_merge <- rbind(g1_tophits[,c("chrpos", "ID", "id", "rsid","LP","ES","SE","REF", "ALT","seqnames","start","end")], 
                      g1_g2tophits[!g1_g2tophits$rsid %in% g1_tophits$rsid, c("chrpos", "ID", "id", "rsid","LP","ES","SE","REF", "ALT","seqnames","start","end")])
    g2_merge <- rbind(g2_tophits[,c("chrpos", "ID", "id", "rsid","LP","ES","SE","REF", "ALT","seqnames","start","end")], 
                      g2_g1tophits[!g2_g1tophits$rsid %in% g2_tophits$rsid, c("chrpos", "ID", "id", "rsid","LP","ES","SE","REF", "ALT","seqnames","start","end")])
    
    write.table(g1_merge, "g1_merge") ## DBEUG
    write.table(g2_merge, "g2_merge")
    
    ##### extract regions for tophit SNPs in either gwas
    regions <- paste0(g1_merge$seqnames, ":", g1_merge$start - 50000, "-", g1_merge$end + 50000)
    
    
    write.table(as.data.frame(regions), "DEBUG_REGIONS_OUTPUT.txt", row.names  =F)
    
    
    regions <- lapply(regions, function(r){
      
      a <- rbind(suppressMessages(vcf_to_tibble(query_gwas(g1_path, chrompos = r))),
                 suppressMessages(vcf_to_tibble(query_gwas(g2_path, chrompos = r))))  %>%
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
          tmp <- x$REF[index]
          x$REF[index] <- x$ALT[index]
          x$ALT[index] <- tmp
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
    
    return(list(g1_raw = g1_merge, g2_raw = g2_merge, regions = regions ))
    
  }, error = function(e) {
    message("Error occurred: ", e$message, " Skipping current trait")
    return(NULL)
  })
  
}

run_fema <- function(betas, ses) {
  
  if (nrow(betas) >=1){
  
  w <- 1 / ses^2
  beta <- rowSums(betas * w) / rowSums(w, na.rm=TRUE)
  se <- sqrt(1 / rowSums(w, na.rm=TRUE))
  z <- abs(beta / se)
  p <- pnorm(z, lower.tail = FALSE)
  nstudy <- apply(betas, 1, \(x) sum(!is.na(x)))
  return(tibble(nstudy, p, z=z))
  
  } else{
    message("too few SNPs for meta analysis")
    
  }
  
}
  

get_region_instruments <- function( instrument_regions, instrument_raw ) {
  
  # slightly hacky fix to name regions
  names(instrument_regions) <- instrument_raw$rsid
  
  
  # remove duplicated ids from instrument_regions
  instrument_regions <- lapply(instrument_regions, \(x) {
    lapply(x, \(y) {
      y[!duplicated(y$rsid), ]
    })
  })
  
  # remove regions that have no data
  rem <- lengths(instrument_regions) == 0
  if (any(rem)) {
    message(sum(rem), " regions have no data")
    instrument_regions <- instrument_regions[!rem]
  }
  
  # for every region, create mats, do scan
  d <- lapply(instrument_regions, \(x) {
    if(is.null(x)) return(NULL)
    if(nrow(x[[1]]) == 0) return(NULL)
    
    rsidintersect <- Reduce(intersect, lapply(x, \(r) r$rsid))
    x <- lapply(x, \(r) {
      r %>% dplyr::filter(rsid %in% rsidintersect) %>% dplyr::filter(!duplicated(rsid)) %>% dplyr::arrange(rsid)
    })
    
    d <- dplyr::select(x[[1]], seqnames, start, rsid, ALT, REF, rsid, id)

    d1 <- run_fema(
        sapply(x, \(y) y$ES), 
        sapply(x, \(y) y$SE)
    )
    
    d <- dplyr::bind_cols(d, d1)
  })
  
  # Keep best SNP from each region
  names(d) <- names(instrument_regions)
  d[sapply(d, is.null)] <- NULL
  d[sapply(d, \(d1) {nrow(d1) == 0})] <- NULL
  dsel <- lapply(d, \(x) {
    subset(x, z == max(x$z, na.rm=TRUE))[1,]
  }) %>% dplyr::bind_rows()
  dsel$region <- names(d)
  # Extract best SNPs from regions for each pop
  inst <- lapply(1:nrow(dsel), \(i) {
    lapply(instrument_regions[[dsel$region[i]]], \(p) {
      
      subset(p, rsid == dsel$rsid[i])
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(id_snp=names(instrument_regions[[dsel$region[i]]]))
  }) %>% dplyr::bind_rows()  %>%
   dplyr::filter(!duplicated(paste(id, rsid)))

  return(inst)
}

  


heterogeneity_calcs <- function(df1, df2, method){
  
  if (method == "raw"){
      
      g1xg2 <- merge(df1[df1$LP > 7.3,], df2, by= "rsid")
      g2xg1 <- merge(df2[df2$LP > 7.3,], df1, by= "rsid")

      g1xg2 <- heterogeneity(g1xg2)
      g2xg1 <- heterogeneity(g2xg1)
        
      out <- as.data.frame(rbind(g1xg2, g2xg1))
      out[,1] <- c(unique(df1$id), unique(df2$id))
      out$method = "raw"
      
  }    
  
  
  if (method == "fema"){
    
    
    
    FEMA_SNPS <- get_region_instruments(df1, df2[[1]])
    
    het_i <- merge(df2[[1]], df2[[2]], by="rsid")
    het_i <- het_i[het_i$rsid %in% FEMA_SNPS$rsid, ]
    
    out <- data.frame(t(heterogeneity(het_i)))
    out[,1] <- unique(df2[[1]]$id)
    out$method = "fema"                       
  }
  
  return(out)
  
}



###########################################################

################## 1. Read in trait list ##################

###########################################################


latAm_gwas <- fread("multi_gwas_across_sources.txt")

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
unique_traits <- unique_traits[5:9] ## DEBUG
for (current_trait in unique_traits)  {
  

  g_l <- latAm_gwas[latAm_gwas$trait == current_trait,]
  g_l <- g_l[order(g_l$sample_size, decreasing = T),]
  
  ids_m <-c(g_l[1,]$id, g_l[2,]$id)
  
  if (file.exists(paste0("vcfs/",ids_m[1], ".vcf.gz")) &
      file.exists(paste0("vcfs/",ids_m[2], ".vcf.gz"))){
    
            p <- pheno_harm(ids_m)
            all_phen <- dplyr::bind_rows(all_phen, p)
            
            instruments <- get_instruments(ids_m)
            
            if (!is.null(instruments)) {
              
              h <- heterogeneity_calcs(instruments$g1_raw, instruments$g2_raw, "raw")
              het <- dplyr::bind_rows(het, h)
              
              h <- heterogeneity_calcs(instruments[[3]], instruments[1:2], "fema")
              het <- dplyr::bind_rows(het, h)
              
              write.table(all_phen, paste0(current_trait, "_pheno_results.txt"), quote = FALSE, row.names = FALSE)
              write.table(het, paste0(current_trait, "_het_results.txt"), quote = FALSE, row.names = FALSE)
              
            } else {
              message("Too few SNPs for comparison.")
            }
            
            
  } else { print("File(s) does not exists")}
}

write.table(all_phen, "pheno_results.txt",quote = F, row.names = F)
write.table(het, "het_results.txt",quote = F, row.names = F)







