library(tidyverse)
setwd("/home/projects/cu_10103/data/INSIGHT_genetics/SMART/INPUT_SMART_05Sep2018")

#import basic phenotype files (including labs) - filter to only include baseline (i.e. visit 0)

SMARTlab <- read_tsv("smart_labs.tsv")
colnames(SMARTlab)
head(SMARTlab)
SMARTlab_base <- filter(SMARTlab, visit == 0)

SMART <- read_tsv("smart.tsv")
SMARTlab_base <- left_join(SMARTlab_base, SMART, by = "pid")

#import PC eigenvectors
setwd("/home/projects/cu_10103/scratch/TET2CGA/eigenstrat_PCA")

PCfile <- read.table("SMART_smartpca_probeQC_recomcutoff.pca.evec", col.names = c("id", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "?"), header = F, skip = 1)

#remove extra id in the id column
library(stringr)

SMARTPCA <- PCfile %>% as_tibble 
SMARTPCA <- separate(SMARTPCA, col = 1, into = c("ID", "ID2"), sep = ":")

#import genotype file and transform to additive format

setwd("/home/projects/cu_10103/scratch/TET2CGA/plink_files")

SNP_screen2tped <- read_delim("/home/projects/cu_10103/scratch/TET2CGA/plink_files/SMART_SNPs_QCed_genericMAF_TET2_IDH12.tped", 
                              "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
SNP_screen2tfam <- read_delim("/home/projects/cu_10103/scratch/TET2CGA/plink_files/SMART_SNPs_QCed_genericMAF_TET2_IDH12.tfam", 
                              "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
SNPmat <- as.matrix(SNP_screen2tped[,-c(1:4)])
rownames(SNPmat) <- as.matrix(SNP_screen2tped[,2])
colnames(SNPmat) <- as.matrix(SNP_screen2tfam[,1])
SNPmat_subset <- SNPmat

SNPmat_subset_ADD <-t(SNPmat_subset)
SNPmat_subset_ADD[SNPmat_subset_ADD == "2 2"] <-0
SNPmat_subset_ADD[SNPmat_subset_ADD == "1 2"] <-1
SNPmat_subset_ADD[SNPmat_subset_ADD == "1 1"] <-2
SNPmat_subset_ADD[SNPmat_subset_ADD == "0 0"] <-NA 
str(SNPmat_subset_ADD)
SNPmat_df <- as.data.frame(SNPmat_subset_ADD) %>% rownames_to_column(., var = "ID")

SMART_geno <- as_tibble(SNPmat_df)

#join tables and create phenogeno tibble

SMARTlab_base$pid <- as.character(SMARTlab_base$pid)

SMARTlab_base <- SMARTlab_base %>%  rename("ID" = "pid")

is.character(SMARTlab_base$pid)



SMART_pheno <- left_join(SMARTlab_base, SMARTPCA, by = "ID")
SMART_genopheno <- left_join(SMART_geno, SMART_pheno, by = "ID")
setwd("/home/projects/cu_10103/scratch/TET2CGA/SMART")
write_tsv(SMART_genopheno, "SMART_genopheno_feb15.tsv")

