# Filter variants based on pre-defined genes from chd_ukbb cohort

library(dplyr)  
library(tidyr)

# Read in data
df <- read.table("~/projects/local/chd_ukbb_archive_032023/chd_ukbb.sample_qced.variant_qced.maf1.vgs.15032023.ht.tsv.bgz", 
                 header = TRUE, 
                 sep = "\t")

#### Filter variants from CHD case-control cohort withing overlapping genes CHD/Cancer ####

# overlap cancer and chd gene list
chd_cancer_overlap <- read.delim('analysis/overlap_chd_cancer_gene_list.tsv') %>% 
  select(gene) %>% 
  pull()

  
# filter dataframe based on gene list
df_filtered <- df %>% 
  filter(SYMBOL %in% chd_cancer_overlap) %>%
  select(-c(csq_group))


# write out filtered dataframe
write.table(df_filtered, 
            file = "analysis/chd_cancer_gene_list_rare_variants_22112023.tsv", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE)


#### Filter variants from CHD case-control cohort withing Cancer genes (germline) ####

# cancer gene list
cancer_genes <- read.delim('./processed_data/genes_list/cosmic/cancer_gene_list_curated.tsv') %>% 
  filter(is_germline == T) %>% 
  select(gene) %>% 
  pull()


# filter dataframe based on gene list
maf = 0.001
df_filtered <- df %>% 
  filter(SYMBOL %in% cancer_genes,
         Consequence != 'synonymous_variant',
         gnomad_genomes_af < maf | is.na(gnomad_genomes_af),
         internal_af < maf | is.na(internal_af)) %>%
  select(-c(csq_group))


# write out filtered dataframe
write.table(df_filtered, 
            file = "analysis/cancer_germline_genes_list_rare_variants_24012024.tsv", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE)



  