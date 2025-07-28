library(gwasglue)
library(ieugwasr)
library(dplyr)
library(readr)
library(data.table)


setwd("WD")
token <- "<TOKEN>"



##Load ancestry and study info data
anc_dat <- read.delim("gwas_cat-anc_info.txt", sep = "\t", stringsAsFactors = FALSE)
traits  <- read.delim("gwas_cat_study_info.txt", sep = "\t", stringsAsFactors = FALSE)

# Merge ancestry + trait info
merged_cat <- merge(traits, anc_dat, by = "STUDY.ACCESSION")

### 2. Pre-filter for relevant studies
filtered_cat <- merged_cat %>%
  filter(`BROAD.ANCESTRAL.CATEGORY` == "Hispanic or Latin American",
         STAGE == "initial",
         !is.na(REPLICATION.SAMPLE.DESCRIPTION),
         `FULL.SUMMARY.STATISTICS` == "yes",
         `NUMBER.OF.INDIVDUALS` >= 1000)

### 3. Collapse to 1 GWAS per PMID + trait (keep max N)
filtered_cat <- filtered_cat %>%
  mutate(sample_size = as.numeric(`NUMBER.OF.INDIVDUALS`),
         trait = MAPPED_TRAIT,
         pmid = PUBMEDID.x,
         id = `STUDY.ACCESSION`,
         author = STUDY,
         ancestry = `BROAD.ANCESTRAL.CATEGORY`) %>%
  group_by(pmid, trait) %>%
  slice_max(order_by = sample_size, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(id, pmid, author, trait, ancestry, sample_size)

### 4. Load and process OpenGWAS metadata
t <- gwasinfo(id = NULL, opengwas_jwt = token)
t <- as.data.frame(t)

Lat_ams <- t %>%
  filter(population == "Hispanic or Latin American",
         sample_size >= 1000,
         trait != "NA") %>%
  select(id, pmid, author, trait, ancestry = population, sample_size) %>%
  distinct(pmid, trait, .keep_all = TRUE)

### 5. Merge both sources
filtered_cat$source <- "EBI"
Lat_ams$source <- "IEU"

combined <- bind_rows(filtered_cat, Lat_ams)


### 6. Identify traits with ≥2 independent studies
traits_with_multiple <- combined %>%
  group_by(trait) %>%
  summarise(n_studies = n()) %>%
  filter(n_studies > 1)

# Filter down to those traits
multi_gwas <- combined %>%
  filter(trait %in% traits_with_multiple$trait)


### 7. Manual filtering
multi_gwas <- multi_gwas %>% 
      filter(!trait %in% c("BMI-adjusted waist-hip ratio", "body height") & ! id %in% c("ebi-a-GCST90095034"))

length(unique(multi_gwas$trait))



### 8. Output
write.table(multi_gwas, "multi_gwas_across_sources.txt", row.names = FALSE, quote = FALSE, sep = "\t")


######################################

## extract matching euro pairs ##

traits <- unique(multi_gwas$trait)

euros <- t %>%
            filter(trait %in% traits & population =="European")

euros_filt <- data.frame()
for (trait in traits){
  
  ts <- euros[euros$trait == trait | tolower(euros$trait) == tolower(trait),]
  
  if (length(ts$id) >= 2){
    ts <- ts[order(as.numeric(ts$sample_size)),]
    ts <- ts[1:2,]
    euros_filt <- dat <- rbind(euros_filt, ts)
  }

}

write.table(euros_filt, "euro_gwas_pairs.txt", row.names = FALSE, quote = FALSE, sep = "\t")

## in total
  # 14 traits
  # 30 GWAS
  # 7 in europeans
