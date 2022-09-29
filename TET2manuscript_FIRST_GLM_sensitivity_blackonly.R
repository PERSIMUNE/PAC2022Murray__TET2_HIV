setwd("/home/projects/cu_10103/scratch/TET2CGA/FIRST")
library(tidyverse)
library(readxl)
library(broom)

FIRST_genopheno <- read_tsv("FIRST_geno_sigonly.tsv")

FIRST_full <- read_tsv("FIRST_genopheno_feb15data.tsv") %>% select("ID", "race", "GENDER", starts_with("PC"),"cd4bl", "rnabl")

FIRST_genopheno <- left_join(FIRST_genopheno, FIRST_full, by = "ID")

##select black only

FIRST_genopheno <- FIRST_genopheno %>% filter(FIRST_genopheno$race == 1)


#set up for GLM

linregtet2 <- FIRST_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs) 
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)
table(linregtet2$genotype)

##check association between groups aka SNPs and variable (VL)

results_rna <- linregtet2 %>% do(fitSNPs = glm(log10(rnabl) ~ genotype + 
                                                 PC1 + PC2 + PC3 +PC4 + 
                                                 GENDER, data = .))
tidyres_rna <- tidy(results_rna, fitSNPs)

tidyres_rna_man <- tidyres_rna %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% 
  rename(rna_pvalue = p.value, rna_estimate = estimate, rna_stderror = std.error)


tidyres_rna_bh <- tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue) %>% 
  mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method="BH"))




###save
write_tsv(tidyres_rna_bh, "FIRST_blackonly_pandqvalue.tsv")

