
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")

# load packages via pacman

pacman::p_load(
  glue,
  KScorrect,
  pscl,
  moments,
  readr,
  DescTools,
  fitdistrplus,
  tidyverse,
  ggplot2,
  ggstatsplot,
  data.table,
  readxl,
  tidyr,
  gridExtra,
  ggpubr,
  MASS,
  rstatix,
  qqman
)

# set work dir

user <- Sys.getenv("USERPROFILE")
setwd(glue(r'({user}\Nextcloud)'))

# remove files from env
# remove all previous plots

rm(list = ls())
graphics.off()


############## Note on file ########################
#                                                 ##
#     Workflow on cancer predisposition genes     ##
#                                                 ##
####################################################

# load CHD-cancer table (raw data)
gene_set_df <-
  read.table(
    glue(r'({user}\Nextcloud\CHD - Cancer\raw_data\variants_table\chd\chd_ukbb.sample_qced.variant_qced.maf1.vgs.15032023.ht.tsv.bgz)'),
    sep = '\t',
    header = T
  )

## load cancer related gene list (from cosmic)
cancer <-
  glue(r'({user}\Nextcloud\CHD - Cancer\raw_data\genes_list\cosmic\Cosmic_Gene_Census_all_18072023.tsv)')
cancer <- read.table(cancer, sep = '\t', header = T)
cancer <- subset(cancer, Germline == 'yes' & Somatic != 'yes')

# load CHD related gene list
# AR = Autosomal Recessive; AD = Autosomal Dominant; XLR = X-linked Recessive; XLD = X-linked Dominant; NA = unknown
chd <-
  glue(r'({user}\Nextcloud\CHD - Cancer\raw_data\genes_list\chd\Gene_list_CHD_reclassification_20201111_FINAL.xlsx)')
chd <- read_excel(chd, sheet = 'CHD gene list_november 2020')

# load CHD_Cancer_overlap list
chd_cancer_overlap <- glue(r'({user}\Nextcloud\CHD - Cancer\analysis\overlap_chd_cancer_gene_list.xlsx)')
chd_cancer_overlap <- read_excel(chd_cancer_overlap, sheet = 'Filtered_Germline')

# split df into het cases
# use allele frequency for filtering
gene_set_df <-
  replace(gene_set_df, is.na(gene_set_df), 0) # replacing NA values with 0s

# filter for syn. variants + allele freq.
cut_off <- 0.001 # Allele frequency cut off

gene_set_df %>%
  filter(n_het_cases > 0 | n_het_controls > 0) %>%
  filter(internal_af < cut_off) %>%
  filter(gnomAD_AF < cut_off) %>% 
  filter(gnomad_genomes_af < cut_off) %>%
  filter(ger_af < cut_off) %>%
  filter(rumc_af < cut_off) %>%
  filter(bonn_af < cut_off) -> df_het

# split df into synonymous, non-synonymous, loss of function (LOF) and missense (all | Constraint) cases based on het. mut.
df_het %>%
  filter(Consequence == "synonymous_variant") -> df_syn

df_het %>%
  filter(Consequence != "synonymous_variant") -> df_nonsyn

df_het[grep("hcLOF,hcLOF_missC", df_het$csq_group), ] -> df_LOF

df_het %>%
  filter(Consequence == "missense_variant") -> df_miss_all

df_het %>%
  filter(Consequence == "missense_variant") %>%
  filter(vep.CADD_PHRED > 20) -> df_miss_C

df <- list("syn" = df_syn,
           "nonsyn" = df_nonsyn,
           "LOF" = df_LOF,
           "all_miss" = df_miss_all,
           "miss_C" = df_miss_C)

t_case <- 3877  # total number of all cases
t_control <- 45083  # total number of all controls

# group by gene name and patient ID -> count
#-------------------------------------------

## data frame selection
## 1 = synonymous, 2 = non-synonymous, 3 = loss of function, 4 = all missense, 5 = missense constraint
d <- 1

df_gene_case <-
  tidyr::separate_rows(data = df[[d]], het_case_ids, sep = "\\|")  # separate chars in cells by | operator
df_gene_case <-
  group_by(df_gene_case, SYMBOL) %>%  # group by SYMBOL
  reframe(het_case_ids = unique(het_case_ids))  # uniquify case IDs

id_count_case <-
  count(df_gene_case, SYMBOL)  # count number of identic SYMBOLS
id_count_case <-
  count(df_gene_case, SYMBOL, het_case_ids != "") # count number of empty cells in case IDs
id_count_case %>%
  filter(!(duplicated(SYMBOL) |
             duplicated(SYMBOL, fromLast = T))) -> nondup_case    # filter df to non duplicated cases
nondup_case$`het_case_ids != ""` <-
  !nondup_case$`het_case_ids != ""`                # invert T to F
nondup_case %>%
  mutate(n = if_else(cumsum(n == 0) > 0, 0, n)) ->  nondup_case                      # set values in column n to 0
id_count_case <-
  rbind(id_count_case, nondup_case)                         # append to original df
id_count_case <-
  id_count_case[!(id_count_case$`het_case_ids != ""` == "FALSE"),]  # remove rows with empty cells
id_count_case <- id_count_case[-2]  # remove second column


df_gene_control <-
  tidyr::separate_rows(data = df[[d]], het_control_ids, sep = "\\|")
df_gene_control <-
  group_by(df_gene_control, SYMBOL) %>%
  reframe(het_control_ids = unique(het_control_ids))
id_count_control <-
  count(df_gene_control, SYMBOL, het_control_ids != "")
id_count_control %>%
  filter(!(duplicated(SYMBOL) |
             duplicated(SYMBOL, fromLast = T))) -> nondup_control
nondup_control$`het_control_ids != ""` <-
  !nondup_control$`het_control_ids != ""`
nondup_control %>%
  mutate(n = if_else(cumsum(n == 0) > 0, 0, n)) ->  nondup_control
id_count_control <- rbind(id_count_control, nondup_control)
id_count_control <-
  id_count_control[!(id_count_control$`het_control_ids != ""` == "FALSE"),]
id_count_control <- id_count_control[-2]


df_gene <- merge.data.frame(
  id_count_case,
  id_count_control,
  by = "SYMBOL"
)

# run fisher test on each gene in relation with patient ID
fisher_test_gene <-
  apply(df_gene, 1,
        function(x) {
          tbl <-
            matrix(
              c(
                as.numeric(x[2]),
                as.numeric(x[3]),
                t_case - as.numeric(x[2]),
                t_control - as.numeric(x[3])
              ),
              nrow = 2,
              ncol = 2,
              byrow = T
            )
          fisher.test(tbl, alternative = "two.sided")
        })

for (i in 1:length(fisher_test_gene)) {
  df_gene[i, "p_fisher"] <- fisher_test_gene[[i]]$p.value
  df_gene[i, "log10"] <- -log10(fisher_test_gene[[i]]$p.value) # add -log10 values for better interpretation of qq-plot
  df_gene[i, "lower"] <- fisher_test_gene[[i]]$conf.int[1]
  df_gene[i, "upper"] <- fisher_test_gene[[i]]$conf.int[2]
  df_gene[i, "or_fisher"] <-
    fisher_test_gene[[i]][["estimate"]][["odds ratio"]]
  df_gene[i, "mean"] <-
    mean(fisher_test_gene[[i]]$conf.int[1], fisher_test_gene[[i]]$conf.int[2])
}

# filter main list for cancer
df_gene %>%
  filter(SYMBOL %in% cancer$Gene.Symbol) -> fisher_cancer

# filter main list for cancer
df_gene %>%
  filter(SYMBOL %in% chd$`genes CHD`) -> fisher_chd
#####

## Plot
#------
n <- names(df[d])

## Plot for all genes in a 2x2 plot env
# par(mfrow = c(2,2))
# qq(df_gene$p_fisher, main = "QQ-Plot of all genes")
# mtext(bquote(paste('het., '*.(n),'., AF'[all]*' < '*.(cut_off))))

## qqplot per gene
## '[all]*','[internal]*','[gnomAD_EX]*','[gnomAD_WG]*','[German]*','[Tueb.]*','[Bonn]*'
par(mfrow = c(1,2))
qq(fisher_cancer$p_fisher, main = "QQ-Plot of cancer predisposition genes")
mtext(bquote(paste('het., '*.(n),'., AF'[all]*' < '*.(cut_off))))

qq(fisher_chd$p_fisher, main = "QQ-Plot of chd genes")
mtext(bquote(paste('het., '*.(n),'., AF'[all]*' < '*.(cut_off))))

#####

## Variant case selection
## 1 = synonymous, 2 = non-synonymous, 3 = loss of function, 4 = missense, 5 = all missense constraint
#------------------------
v <- 3
#####

## Cancer
#--------
df[[1]] %>%
  filter(SYMBOL %in% cancer$Gene.Symbol) -> cancer_syn
df[[2]] %>%
  filter(SYMBOL %in% cancer$Gene.Symbol) -> cancer_nonsyn
df[[3]] %>%
  filter(SYMBOL %in% cancer$Gene.Symbol) -> cancer_LOF
df[[4]] %>%
  filter(SYMBOL %in% cancer$Gene.Symbol) -> cancer_miss_all
df[[5]] %>%
  filter(SYMBOL %in% cancer$Gene.Symbol) -> cancer_miss_C

var_cancer <- list("cancer_syn" = cancer_syn,
                   "cancer_nonsyn" = cancer_nonsyn,
                   "cancer_LOF" = cancer_LOF,
                   "cancer_miss-all" = cancer_miss_all,
                   "cancer_miss-C" = cancer_miss_C)

# break down to unique ids
# sum up het. cases and het. controls not group by gene
# make fisher exact test

cancer_case <-
   tidyr::separate_rows(data = var_cancer[[v]], het_case_ids, sep = "\\|")  %>% # separate chars in cells by | operator
   filter(startsWith(x = het_case_ids, "BIID"))
cancer_case <-
  group_by(cancer_case, het_case_ids) %>%
  reframe(het_case_ids = unique(het_case_ids)) # uniquify case IDs
n_het_case_ID_cancer <- length(cancer_case$het_case_ids)

cancer_control <-
   tidyr::separate_rows(data = var_cancer[[v]], het_control_ids, sep = "\\|") %>%
   filter(startsWith(x = het_control_ids, "BIID"))
cancer_control <-
  group_by(cancer_control, het_control_ids) %>%
  reframe(het_control_ids = unique(het_control_ids)) # uniquify case IDs
n_het_control_ID_cancer <- length(cancer_control$het_control_ids)

m_case_control_cancer <-
  matrix(
    c(
      n_het_case_ID_cancer,
      n_het_control_ID_cancer,
      t_case - n_het_case_ID_cancer,
      t_control - n_het_control_ID_cancer
    ),
    nrow = 2
  )
  
het_fisher_cancer <- 
  fisher.test(m_case_control_cancer, alternative = "two.sided")
#####

## CHD
#-----
df[[1]] %>%
  filter(SYMBOL %in% chd$`genes CHD`) -> chd_syn
df[[2]] %>%
  filter(SYMBOL %in% chd$`genes CHD`) -> chd_nonsyn
df[[3]] %>%
  filter(SYMBOL %in% chd$`genes CHD`) -> chd_LOF
df[[4]] %>%
  filter(SYMBOL %in% chd$`genes CHD`) -> chd_miss_all
df[[5]] %>%
  filter(SYMBOL %in% chd$`genes CHD`) -> chd_miss_C

var_chd <- list("chd_syn" = chd_syn,
                "chd_nonsyn" = chd_nonsyn,
                "chd_LOF" = chd_LOF,
                "chd_miss-all" = chd_miss_all,
                "chd_miss-C" = chd_miss_C)

chd_case <-
   tidyr::separate_rows(data = var_chd[[v]], het_case_ids, sep = "\\|") %>% # separate chars in cells by | operator
   filter(startsWith(x = het_case_ids, "BIID"))
chd_case <-
  group_by(chd_case, het_case_ids) %>%
  reframe(het_case_ids = unique(het_case_ids)) # uniquify case IDs
n_het_case_ID_chd <- length(chd_case$het_case_ids)

chd_control <-
   tidyr::separate_rows(data = var_chd[[v]], het_control_ids, sep = "\\|") %>%
   filter(startsWith(x = het_control_ids, "BIID"))
chd_control <-
  group_by(chd_control, het_control_ids) %>%
  reframe(het_control_ids = unique(het_control_ids)) # uniquify case IDs
n_het_control_ID_chd <- length(chd_control$het_control_ids)


m_case_control_chd <-
  matrix(
    c(
      n_het_case_ID_chd,
      n_het_control_ID_chd,
      t_case - n_het_case_ID_chd,
      t_control - n_het_control_ID_chd
    ),
    nrow = 2
  )

het_fisher_chd <- 
  fisher.test(m_case_control_chd, alternative = "two.sided")
#####

## Overlap
#---------
df[[1]] %>%
  filter(SYMBOL %in% chd_cancer_overlap$gene) -> overlap_syn
df[[2]] %>%
  filter(SYMBOL %in% chd_cancer_overlap$gene) -> overlap_nonsyn
df[[3]] %>%
  filter(SYMBOL %in% chd_cancer_overlap$gene) -> overlap_LOF
df[[4]] %>%
  filter(SYMBOL %in% chd_cancer_overlap$gene) -> overlap_miss_all
df[[5]] %>%
  filter(SYMBOL %in% chd_cancer_overlap$gene) -> overlap_miss_C

var_overlap <- list("overlap_syn" = overlap_syn,
                    "overlap_nonsyn" = overlap_nonsyn,
                    "overlap_LOF" = overlap_LOF,
                    "overlap_miss-all" = overlap_miss_all,
                    "overlap_miss-C" = overlap_miss_C)

# break down to unique ids
# sum up het. cases and het. controls not group by gene
# make fisher exact test

overlap_case <-
   tidyr::separate_rows(data = var_overlap[[v]], het_case_ids, sep = "\\|") %>% # separate chars in cells by | operator
   filter(startsWith(x = het_case_ids, "BIID"))
overlap_case <-
  group_by(overlap_case, het_case_ids) %>%
  reframe(het_case_ids = unique(het_case_ids)) # uniquify case IDs
n_het_case_ID_overlap <- length(overlap_case$het_case_ids)

overlap_control <-
   tidyr::separate_rows(data = var_overlap[[v]], het_control_ids, sep = "\\|") %>%
   filter(startsWith(x = het_control_ids, "BIID"))
overlap_control <-
  group_by(overlap_control, het_control_ids) %>%
  reframe(het_control_ids = unique(het_control_ids)) # uniquify case IDs
n_het_control_ID_overlap <- length(overlap_control$het_control_ids)

m_case_control_overlap <-
  matrix(
    c(
      n_het_case_ID_overlap,
      n_het_control_ID_overlap,
      t_case - n_het_case_ID_overlap,
      t_control - n_het_control_ID_overlap
    ),
    nrow = 2
  )

het_fisher_overlap <- 
  fisher.test(m_case_control_overlap, alternative = "two.sided")
#####

## fisher df
#--------------


Chd <- c(het_fisher_chd$p.value,
         het_fisher_chd$conf.int[1],
         het_fisher_chd$conf.int[2],
         het_fisher_chd[["estimate"]][["odds ratio"]])
Cancer <- c(het_fisher_cancer$p.value,
            het_fisher_cancer$conf.int[1],
            het_fisher_cancer$conf.int[2],
            het_fisher_cancer[["estimate"]][["odds ratio"]])
Overlap <- c(het_fisher_overlap$p.value,
             het_fisher_overlap$conf.int[1],
             het_fisher_overlap$conf.int[2],
             het_fisher_overlap[["estimate"]][["odds ratio"]])
df_chd_cancer <- as.data.frame(rbind(Chd,Cancer,Overlap))
df_chd_cancer$Groups <- c("CHD", "Cancer", "Overlap")
colnames(df_chd_cancer) <- c("p_value", "lower", "upper", "OR", "Groups")
#####

## plot odds ratio chd,cancer,overlap
#------------------------------------
# split list names -> use for subtitle
s <- strsplit(names(var_cancer[v]), "_", perl = T)[[1]][[2]] # change number according to var case

f_plot <-
  ggplot(df_chd_cancer, aes(y = Groups)) +
  theme_classic() +
  geom_point(aes(x = OR), shape = 19, size = 2) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = 'red') +
  labs(x = "OR", y = "",
       title = "Odds ratio plot of the following gene-set:",
       subtitle = bquote(paste('Cancer predisposition, CHD, Overlap (CHD-Cancer), het., '*.(s),'., AF < '*.(cut_off)))) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 7)
  )
f_plot

#####













