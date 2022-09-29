library(tidyverse)
setwd("/home/projects/cu_10103/data/INSIGHT_genetics/START/START_27102017")

# Import START phenotype data
STARTgenomics <- read_tsv("startgenomics.tsv")

START_basic <- select(STARTgenomics, "ID", "GENDER", "age_00", starts_with("wbc_c"), "cgroup", "nowsmoke_00", "exsmoke_00", "race", "rtreat", "RNA_00", "startedrna", "cd4nadir_00", "CD4_00", "cd4cd8_00", "WBC_00", "LYMPHO_00", "CD8_00", "cgroup", "RNA_01", "RNA_04", "RNA_08", "RNA_12")


setwd("/home/projects/cu_10103/data/INSIGHT_genetics/START/START_24042018")
STARTbiomarker <- read_tsv("biomarkers.tsv")
STARTbiomarker_base <- filter(STARTbiomarker, VISIT == 0)
START_pheno <- left_join(START_basic, STARTbiomarker_base, by = "ID")



# import PC eigenvectors
setwd("/home/projects/cu_10103/scratch/TET2CGA/eigenstrat_PCA")

PCfile <- read.table("START_smartpca_probeQC_recomcutoff.pca.evec", col.names = c("id", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "?"), header = F, skip = 1)
# merge pheno and PC file


# remove extra id in the id column
library(stringr)

STARTPCA <- PCfile %>% as_tibble()
STARTPCA <- separate(STARTPCA, col = 1, into = c("ID", "ID2"), sep = ":")

class(STARTPCA$ID)
START_pheno$ID <- as.character(START_pheno$ID)
class(START_pheno$ID)
as.character(START_pheno$ID)
START_pheno_PC <- left_join(START_pheno, STARTPCA, by = "ID")

# import genotype file and transform to additive format

setwd("/home/projects/cu_10103/scratch/TET2CGA/plink_files")

SNP_screen2tped <- read_delim("/home/projects/cu_10103/scratch/TET2CGA/plink_files/START_SNPs_QCed_genericMAF_TET2_IDH12.tped",
  "\t",
  escape_double = FALSE, col_names = FALSE,
  trim_ws = TRUE
)
SNP_screen2tfam <- read_delim("/home/projects/cu_10103/scratch/TET2CGA/plink_files/START_SNPs_QCed_genericMAF_TET2_IDH12.tfam",
  "\t",
  escape_double = FALSE, col_names = FALSE,
  trim_ws = TRUE
)
SNPmat <- as.matrix(SNP_screen2tped[, -c(1:4)])
rownames(SNPmat) <- as.matrix(SNP_screen2tped[, 2])
colnames(SNPmat) <- as.matrix(SNP_screen2tfam[, 1])
SNPmat_subset <- SNPmat

SNPmat_subset_ADD <- t(SNPmat_subset)
SNPmat_subset_ADD[SNPmat_subset_ADD == "2 2"] <- 0
SNPmat_subset_ADD[SNPmat_subset_ADD == "1 2"] <- 1
SNPmat_subset_ADD[SNPmat_subset_ADD == "1 1"] <- 2
SNPmat_subset_ADD[SNPmat_subset_ADD == "0 0"] <- NA
str(SNPmat_subset_ADD)
SNPmat_df <- as.data.frame(SNPmat_subset_ADD) %>% rownames_to_column(., var = "ID")

START_geno <- as_tibble(SNPmat_df)

# join tables and create phenogeno tibble

START_pheno$ID <- as.character(START_pheno$ID)
is.character(START_pheno$ID)
is.character(STARTPCA$ID)

START_pheno <- left_join(START_pheno, STARTPCA, by = "ID")
START_genopheno <- left_join(START_geno, START_pheno, by = "ID")
setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
write_tsv(START_genopheno, "START_genopheno_15febdata.tsv")
str(START_genopheno$CD4_00)


## subset main geno df with significnat SNPs from the START and FIRST associations
# select SNPs with q value < 0.05 from START and FIRST associations

START_geno_sigonly <-
  START_geno %>% select(
    "ID",
    "AX-14431677",
    "AX-175346011",
    "AX-14431796",
    "AX-34492519",
    "AX-151531994",
    "AX-14431777",
    "AX-14431755",
    "AX-14431884",
    "AX-14431750",
    "AX-14431754",
    "AX-14431783",
    "AX-41321305",
    "AX-14431833",
    "AX-34492557",
    "AX-175364028",
    "AX-38155539",
    "AX-14431819",
    "AX-14431848",
    "AX-14431808",
    "AX-14431631",
    "AX-14431817",
    "AX-175367051",
    "AX-14431804",
    "AX-14431806",
    "AX-14431861",
    "AX-14431887",
    "AX-175346090",
    "AX-14431860",
    "AX-14431883",
    "AX-14431884",
    "AX-14431860",
    "AX-14431848",
    "AX-14431804",
    "AX-14431806",
    "AX-14431783",
    "AX-14431872",
    "AX-14431808",
    "AX-14431819",
    "AX-175346090",
    "AX-14431754",
    "AX-151531994",
    "AX-34492557",
    "AX-14431847",
    #    "AX-14431803",
    #    "AX-14431813",
    "AX-33487697",
    #   "AX-34492391",
    "AX-175345974"
  )




setwd("/home/projects/cu_10103/scratch/TET2CGA/Descriptive tables")
write_tsv(START_geno_sigonly, "START_geno_sigonly.tsv")




###merge with cancer events for COX models
###Note, the cancer events are not relevant for the manuscript titled Association between TET2 genetic variation and viral load in people living with HIV
###But the genopheno file created is referenced in the sensittivity analysis code and it is therefore important to show how this was created

setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
cancer <- read_tsv("START_mds_cancer_events.tsv") 
cancer$ID <- as.character(cancer$ID)

cancer_pheno <- full_join(cancer, START_pheno_PC, by = "ID")

count(cancer_pheno, factor(cancer_pheno$MDS_ever))

###merge with significant SNPs from START and FIRST

setwd("/home/projects/cu_10103/scratch/TET2CGA/Descriptive tables")
sigonly <- read_tsv("START_geno_sigonly.tsv")


#merge sig genopheno
sigonly$ID <- as.character(sigonly$ID)
START_sigsnps_cancer <- full_join(cancer_pheno, sigonly, by = "ID")

setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
write_tsv(START_sigsnps_cancer, "START_sigsnps_cancer_genopheno.tsv")
