# Parse Cosmic Cancer Mutation File

library(dplyr)
library(tidyr)


# import variant file (tsv)
cosmic_variants <- read.delim(file = "raw_data/variants_table/cosmic/CancerMutationCensus_AllData_v98_GRCh38.tsv.gz",
                              sep = "\t")


# Get list of cancer genes
cancer_genes <- read.delim(file = 'processed_data/genes_list/cancer/cancer_gene_list_curated.tsv') %>%
  pull(gene)


# Filter to cancer gene and rare variants in GNOMAD
af_threshold = 0.01

cosmic_variants %>%
  filter(GNOMAD_GENOMES_AF < af_threshold | is.na(GNOMAD_GENOMES_AF),
         GNOMAD_EXOMES_AF < af_threshold | is.na(GNOMAD_EXOMES_AF),
         GENE_NAME %in% cancer_genes) -> df_rare


# write table
write.table(df_rare,
            file = "./processed_data/variants_table/cosmic/CancerMutationCensus_rare_variants_v98_GRCh38.tsv",
            sep = "\t",
            row.names = F,
            quote = F)
