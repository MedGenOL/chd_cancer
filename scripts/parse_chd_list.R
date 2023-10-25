# Script to parse the list of CHD genes

# import libraries
library(dplyr)
library(tidyr)
library(openxlsx)


# read excel file with CHD genes
chd_genes <- read.xlsx("./data/genes_list/chd/Gene_list_CHD_reclassification_20201111_FINAL.xlsx", sheet = 1, 
                        colNames = TRUE, 
                        detectDates = FALSE, 
                        skipEmptyRows = TRUE)

# rename/selec cols
chd_genes %>%
  dplyr::rename(gene = genes.CHD,
                tier_chd = new.ranking.JB) %>%
  select(gene, tier_chd) -> chd_genes

# re-code tier info 
chd_genes %>%
  mutate(tier_chd = case_when(startsWith(tier_chd, "CHD_1") ~ 1,
                              startsWith(tier_chd, "CHD_2") ~ 2,
                              startsWith(tier_chd, "CHD_3") ~ 3,
                              startsWith(tier_chd, "CHD_4") ~ 4,
                              T ~ -1)
         ) -> chd_genes

write.table(chd_genes,
            file = "./data/genes_list/chd/chd_gene_list_curated.tsv",
            sep = "\t",
            row.names = F,
            quote = F)
