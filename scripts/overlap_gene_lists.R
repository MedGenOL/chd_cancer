# Merge list of CHD and cancer genes

library(dplyr)


chd_genes <- read.delim(file = "./data/genes_list/chd/chd_gene_list_curated.tsv")
cancer_genes <- read.delim(file = "./data/genes_list/cancer/cancer_gene_list_curated.tsv")


df_merged <- merge.data.frame(chd_genes,
                              cancer_genes,
                              by = 'gene',
                              all = F)

write.table(df_merged,
            file = "./data/genes_list/overlap/overlap_chd_cancer_gene_list.tsv",
            sep = "\t",
            row.names = F,
            quote = F)
