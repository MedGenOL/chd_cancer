# Filter variants based on pre-defined genes from chd_ukbb cohort

library(dplyr)  
library(tidyr)

# Read in data
df <- read.table("~/projects/local/chd_ukbb_archive_032023/chd_ukbb.sample_qced.variant_qced.maf1.vgs.15032023.ht.tsv.bgz", 
                 header = TRUE, 
                 sep = "\t")

# Read in gene list
genes <- read.delim('analysis/overlap_chd_cancer_gene_list.tsv') %>% 
  select(gene) %>% 
  pull()

  
# filter dataframe based on gene list
df_filtered <- df %>% 
  filter(SYMBOL %in% genes) %>%
  select(-c(csq_group))


# write out filtered dataframe
write.table(df_filtered, 
            file = "analysis/chd_cancer_gene_list_rare_variants_22112023.tsv", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE)



  