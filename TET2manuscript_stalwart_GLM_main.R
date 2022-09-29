setwd("/home/projects/cu_10103/scratch/TET2CGA/STALWART")
library(tidyverse)
library(readxl)
library(broom)

STALWART_genopheno <- read_tsv("STALWART_genopheno_feb15data.tsv")

linregtet2 <- STALWART_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs) 
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)
table(linregtet2$genotype)

##check association between groups aka SNPs and variable (VL)

results_rna <- linregtet2 %>% do(fitSNPs = glm(log10(RNAbl) ~ genotype + PC1 + PC2 + PC3 +PC4 + Gender, data = .))
tidyres_rna <- tidy(results_rna, fitSNPs)

tidyres_rna_man <- tidyres_rna %>%
  filter(term == "genotype") %>% select("SNPs", "p.value") %>% rename(rna_pvalue = p.value)

tidyres_rna_bh <- tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue) %>% 
  mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method="BH"))


##check association between groups aka SNPs and variable (CD4)
### Plus 1 to log transformed data to account for 0's in data
class(STALWART_genopheno$cd4bl)

results_cd4 <- linregtet2 %>% do(fitSNPs = glm(Cd4Bl ~ genotype + PC1 + PC2 + PC3 +PC4 + Gender, data = .))
tidyres_cd4 <- tidy(results_cd4, fitSNPs)

tidyres_cd4_man <- tidyres_cd4 %>%
  filter(term == "genotype") %>% select("SNPs", "p.value") %>% rename(cd4_pvalue = p.value)

tidyres_cd4_bh <- tidyres_cd4_man %>% ungroup(tidyres_cd4_man$cd4_palue) %>% 
  mutate(cd4_qvalue = p.adjust(tidyres_cd4_man$cd4_pvalue, method="BH"))



#combine all pvalues into a signle tibble
getwd()
STALWART_pandqvalues <- left_join(tidyres_cd4_bh, tidyres_rna_bh, by = "SNPs") 

write_tsv(STALWART_pandqvalues, "STALWART_pandqvalues")

