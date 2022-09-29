setwd("/home/projects/cu_10103/scratch/TET2CGA/FIRST")
library(tidyverse)
library(readxl)
library(broom)

FIRST_genopheno <- read_tsv("FIRST_genopheno_feb15data.tsv")
### Create format for multiple linear regression

## selecting and grouping SNPs
linregtet2 <- FIRST_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs) 
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)

## formatting data types
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)

## visually inspect result
table(linregtet2$genotype)

##check association between groups aka SNPs and variable (VL)

results_rna <- linregtet2 %>% do(fitSNPs = glm(log10(rnabl) ~ genotype + PC1 + PC2 + PC3 +PC4 + GENDER, data = .))
tidyres_rna <- tidy(results_rna, fitSNPs)

tidyres_rna_man <- tidyres_rna %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% 
  rename(rna_pvalue = p.value, rna_estimate = estimate, rna_stderror = std.error)


tidyres_rna_bh <- tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue) %>% 
  mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method="BH"))



##check association between groups aka SNPs and variable (CD4)
### Plus 1 to log transformed data to account for 0's in data
class(FIRST_genopheno$cd4bl)

results_cd4 <- linregtet2 %>% do(fitSNPs = glm(cd4bl ~ genotype + PC1 + PC2 + PC3 +PC4 + GENDER, data = .))
tidyres_cd4 <- tidy(results_cd4, fitSNPs)

tidyres_cd4_man <- tidyres_cd4 %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% 
  rename(cd4_pvalue = p.value, cd4_estimate = estimate, cd4_stderror = std.error)


tidyres_cd4_bh <- tidyres_cd4_man %>% ungroup(tidyres_cd4_man$cd4_palue) %>% 
  mutate(cd4_qvalue = p.adjust(tidyres_cd4_man$cd4_pvalue, method="BH"))



#Single regression - check a single regression to confirm association

# 
# summary(glm(log10(rnabl) ~ factor("AX-107670385") + GENDER + PC1 + PC2 + PC3 + PC4, data = FIRST_genopheno)
# 
# tidy(glm(log10(rnabl) ~ factor("AX-107670385") + race + PC1 + PC2 + PC3 + PC4, data=FIRST_genopheno))

#plot 

linregtet2 %>% filter(SNPs == "AX-34492519") %>% ggplot(aes(x=factor(genotype), y=(log10(rnabl)))) + geom_boxplot() + geom_point() + theme_minimal()

#combine all pvalues into a signle tibble
getwd()
FIRST_pandqvalues <- left_join(tidyres_cd4_bh, tidyres_rna_bh, by = "SNPs") 

write_tsv(FIRST_pandqvalues, "FIRST_pandqvalues_feb15")
                                                         