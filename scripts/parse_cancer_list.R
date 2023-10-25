# Script to parse the list of cancer genes

# import libraries
library(dplyr)
library(tidyr)


# import gene list
cancer_genes <- read.delim(file = './data/genes_list/cancer/Cosmic_Gene_Census_all_18072023.tsv',
                           sep = '\t')


# rename/select cols
cancer_genes %>%
  dplyr::select(gene = Gene.Symbol,
                tier_cancer = Tier,
                is_germline = Germline) -> cancer_genes

# re-code germline info 
cancer_genes %>%
  mutate(is_germline = if_else(is_germline == 'yes', T, F)
  ) -> cancer_genes

write.table(cancer_genes,
            file = "./data/genes_list/cancer/cancer_gene_list_curated.tsv",
            sep = "\t",
            row.names = F,
            quote = F)