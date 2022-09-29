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
setwd("/home/projects/cu_10103/scratch/TET2CGA/START")

start <- read_tsv("START_sigsnps_cancer_genopheno.tsv")
as_tibble(start)

start_subset <- dplyr::select(start, c("ID", starts_with("AX"))) %>% column_to_rownames("ID")

start_mat <- as.matrix(start_subset) ###convert to matrix


  # convert to genotype format
  

str(start_subset)


start_allele <- makeGenotypes(start_mat, c(1:33), method = as.genotype.allele.count) 

#reorder columns based on estimate in the START glm - taken from the sigSNPS from START_FIRST_smartfollowup file 
##remove the SNPs not present in START before attempting to reorder
###if changing this and including genetic position, make sure to change the geneposition vector as well


col_order <- c(
  "AX-34492519",
  "AX-14431796",
  "AX-14431677",
  "AX-175346011",
  "AX-14431777",
  "AX-14431750",
  "AX-41321305",
  "AX-175364028",
  "AX-14431833",
  "AX-175367051",
  "AX-14431887",
  "AX-14431631",
  "AX-38155539",
  "AX-33487697",
  "AX-175345974",
  "AX-14431883",
  "AX-14431847",
  "AX-14431817",
  "AX-14431872",
  "AX-34492557",
  "AX-14431860",
  "AX-14431754",
  "AX-14431884",
  "AX-175346090",
  "AX-151531994",
  "AX-14431861",
  "AX-14431806",
  "AX-14431804",
  "AX-14431848",
  "AX-14431808",
  "AX-14431819",
  "AX-14431783",
  "AX-14431755"
)



start_allele_reorder <- start_allele[, col_order]

colnames(start_allele_reorder)

#rename AX IDS with rs-IDS
col_order_df <- as_tibble(col_order) %>% rename("SNPs" = "value")
sigsnpnames <- left_join(col_order_df, sigsnps, by = "SNPs") %>% dplyr::select(`dbSNP RS ID`)

sigsnpnames_vec <- sigsnpnames[[1]]

start_allele_reorder_rsID <- start_allele_reorder %>% rename_at(vars(col_order), ~ sigsnpnames_vec)


# g <- as.genotype.allele.count(start_geno, alleles=c("C","T"))

is.genotype(start_mat
            )
is.genotype(start_allele)



##Perfrom LD

startld <- LD(start_allele_reorder,
               stats=c("D.prime", "R.squared"))


sum <- summary(startld, digits = 3,
           which = c("r"),
           show.all = TRUE)


print(sum)
print(startld, digits = getOption("digits"))

###heatmap


snpnames <- colnames(start_allele_reorder)

#get genetic location in same order as snpnames

start_geneposition <- c(106173177,
                        106151781,
                        106095219,
                        106122657,
                        106141445,
                        106132312,
                        106110796,
                        106162216,
                        106173540,
                        106184236,
                        106200234,
                        106064349,
                        106065102,
                        209128892,
                        106082473,
                        106197750,
                        106182132,
                        106166139,
                        106194419,
                        106190732,
                        106188548,
                        106133471,
                        106197936,
                        106201653,
                        106142080,
                        106189156,
                        106160161,
                        106159386,
                        106182446,
                        106160495,
                        106168623,
                        106145468,
                        106133750
                        )

##perform LDheatmap

dheatmap <- LDheatmap(start_allele_reorder_rsID, genetic.distances=NULL, distances="physical",
          LDmeasure="r", title="START Pairwise LD", add.map=FALSE, add.key=TRUE,
          geneMapLocation=0.15, geneMapLabelX=NULL, geneMapLabelY=NULL,
          SNP.name=sigsnpnames_vec, color=heat.colors(20), newpage=TRUE,
          name="ldheatmap", vp.name=NULL, pop=FALSE, flip=FALSE, text=FALSE)

#extract the triangle matrix with LD (r2) values
dmatrix <- dheatmap$LDmatrix
# then convert the ld matrix to a data frame
mat <- as_tibble(dmatrix) %>% round(3) 

#add snpnames as a column
mat_ <- cbind(snpnames, mat)
getwd()
write_tsv(mat_, "start_sigsnpsLD.tsv")
