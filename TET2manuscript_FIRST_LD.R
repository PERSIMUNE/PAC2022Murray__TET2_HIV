library(tidyverse)
options(pillar.sigfig = 3)
library(dplyr)
library(genetics)
library(broom)
library(ggrepel)
library(ggplot2)

library(LDheatmap)

#load rsIDs so I can merge them into figures more easily
setwd("/home/projects/cu_10103/scratch/TET2CGA/data")
sigsnps <- read_tsv("significantSNPs_fromSTARTandFIRST.tsv")

sigsnps <- sigsnps %>% dplyr::select(`Probe Set ID`, `dbSNP RS ID`) %>% dplyr::rename("SNPs" = "Probe Set ID")

#select SNPs of interest for LD analysis
setwd("/home/projects/cu_10103/scratch/TET2CGA/Descriptive tables")
list.files()
first <- read_tsv("FIRST_geno_sigonly.tsv")
as_tibble(first)

first_subset <- dplyr::select(first, c("ID", starts_with("AX"))) %>% column_to_rownames("ID")

first_mat <- as.matrix(first_subset) ###convert to matrix


  # convert to genotype format
  

str(first_subset)
colnames(first_mat)

first_allele <- makeGenotypes(first_mat, c(1:35), method = as.genotype.allele.count) 

#reorder columns based on estimate in the first glm - taken from the sigSNPS from START_FIRST_smartfollowup file 
##remove the SNPs not present in first before attempting to reorder
###if changing this and including genetic position, make sure to change the geneposition vector as well


col_order <- c(
  "AX-38155539",
  "AX-14431631",
  "AX-175367051",
  "AX-14431677",
  "AX-33487697",
  "AX-34492519",
  "AX-175364028",
  "AX-14431750",
  "AX-175346011",
  "AX-14431796",
  "AX-14431887",
  "AX-41321305",
  "AX-14431777",
  "AX-14431833",
  "AX-175345974",
  "AX-34492391",
  "AX-14431847",
  "AX-14431754",
  "AX-151531994",
  "AX-34492557",
  "AX-14431883",
  "AX-14431819",
  "AX-175346090",
  "AX-14431848",
  "AX-14431808",
  "AX-14431803",
  "AX-14431813",
  "AX-14431872",
  "AX-14431783",
  "AX-14431884",
  "AX-14431804",
  "AX-14431806",
  "AX-14431817",
  "AX-14431861",
  "AX-14431860"
  )



first_allele_reorder <- first_allele[, col_order]

colnames(first_allele_reorder)

#rename AX IDS with rs-IDS
col_order_df <- as_tibble(col_order) %>% rename("SNPs" = "value")
sigsnpnames <- left_join(col_order_df, sigsnps, by = "SNPs") %>% dplyr::select(`dbSNP RS ID`)

sigsnpnames_vec <- sigsnpnames[[1]]

first_allele_reorder_rsID <- first_allele_reorder %>% rename_at(vars(col_order), ~ sigsnpnames_vec)



# g <- as.genotype.allele.count(first_geno, alleles=c("C","T"))

is.genotype(first_mat
            )
is.genotype(first_allele)



##Perfrom LD

firstld <- LD(first_allele_reorder,
               stats=c("D.prime", "R.squared"))


sum <- summary(firstld, digits = 3,
           which = c("r"),
           show.all = TRUE)


print(sum)
print(firstld, digits = getOption("digits"))

###heatmap


snpnames <- colnames(first_allele_reorder)

#get genetic location in same order as snpnames

first_geneposition <- c(106065102,
                        106064349,
                        106184236,
                        106095219,
                        209128892,
                        106173177,
                        106162216,
                        106132312,
                        106122657,
                        106151781,
                        106200234,
                        106110796,
                        106141445,
                        106173540,
                        106082473,
                        106110032,
                        106182132,
                        106133471,
                        106142080,
                        106190732,
                        106197750,
                        106168623,
                        106201653,
                        106182446,
                        106160495,
                        106158738,
                        106164723,
                        106194419,
                        106145468,
                        106197936,
                        106159386,
                        106160161,
                        106166139,
                        106189156,
                        106188548
                                              )

##perform LDheatmap

dheatmap <- LDheatmap(first_allele_reorder_rsID, genetic.distances=NULL, distances="physical",
          LDmeasure="r", title="FIRST Pairwise LD", add.map=FALSE, add.key=TRUE,
          geneMapLocation=0.15, geneMapLabelX=NULL, geneMapLabelY=NULL,
          SNP.name=sigsnpnames_vec, color=heat.colors(20), newpage=TRUE,
          name="ldheatmap", vp.name=NULL, pop=FALSE, flip=FALSE, text=FALSE)



