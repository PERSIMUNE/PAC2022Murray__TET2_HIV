setwd("/home/projects/cu_10103/scratch/TET2CGA/ESPRIT")
library(tidyverse)
library(readxl)
library(broom)

ESPRIT_genopheno <- read_tsv("ESPRIT_genopheno_feb15data.tsv")

linregtet2 <- ESPRIT_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs) 
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)
table(linregtet2$genotype)

##check association between groups aka SNPs and variable (VL)
### Remove as ESPRIT patients are not ART naive

#results_rna <- linregtet2 %>% do(fitSNPs = glm(log10(RNAbl) ~ genotype + PC1 + PC2 + PC3 +PC4 + Gender, data = .))
#tidyres_rna <- tidy(results_rna, fitSNPs)

#tidyres_rna_man <- tidyres_rna %>%
 # filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% rename(rna_pvalue = p.value)

#tidyres_rna_bh <- tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue) %>% 
 # mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method="BH"))


##check association between groups aka SNPs and variable (IL6)

results_il6 <- linregtet2 %>% do(fitSNPs = glm(log2(IL6bl) ~ genotype + PC1 + PC2 + PC3 +PC4 + Gender, data = .))
tidyres_il6 <- tidy(results_il6, fitSNPs)

tidyres_il6_man <- tidyres_il6 %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% rename(il6_pvalue = p.value)

tidyres_il6_bh <- tidyres_il6_man %>% ungroup(tidyres_il6_man$il6_palue) %>% 
  mutate(il6_qvalue = p.adjust(tidyres_il6_man$il6_pvalue, method="BH"))


##check association between groups aka SNPs and variable (CD4)

results_cd4 <- linregtet2 %>% do(fitSNPs = glm(CD4bl ~ genotype + PC1 + PC2 + PC3 +PC4 + Gender, data = .))
tidyres_cd4 <- tidy(results_cd4, fitSNPs)

tidyres_cd4_man <- tidyres_cd4 %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% rename(cd4_pvalue = p.value)

tidyres_cd4_bh <- tidyres_cd4_man %>% ungroup(tidyres_cd4_man$cd4_palue) %>% 
  mutate(cd4_qvalue = p.adjust(tidyres_cd4_man$cd4_pvalue, method="BH"))


##check association between groups aka SNPs and variable (crp)

results_crp <- linregtet2 %>% do(fitSNPs = glm(log2(CRPbl) ~ genotype + PC1 + PC2 + PC3 +PC4 + Gender, data = .))
tidyres_crp <- tidy(results_crp, fitSNPs)

tidyres_crp_man <- tidyres_crp %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% rename(crp_pvalue = p.value)

tidyres_crp_bh <- tidyres_crp_man %>% ungroup(tidyres_crp_man$crp_palue) %>% 
  mutate(crp_qvalue = p.adjust(tidyres_crp_man$crp_pvalue, method="BH"))


##check association between groups aka SNPs and variable (ddimer)

results_ddimer <- linregtet2 %>% do(fitSNPs = glm(log2(ddimerBL) ~ genotype + PC1 + PC2 + PC3 +PC4 + Gender, data = .))
tidyres_ddimer <- tidy(results_ddimer, fitSNPs)

tidyres_ddimer_man <- tidyres_ddimer %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% rename(ddimer_pvalue = p.value)

tidyres_ddimer_bh <- tidyres_ddimer_man %>% ungroup(tidyres_ddimer_man$ddimer_palue) %>% 
  mutate(ddimer_qvalue = p.adjust(tidyres_ddimer_man$ddimer_pvalue, method="BH"))




#plot 

linregtet2 %>% filter(SNPs == "AX-14431796") %>% ggplot(aes(x=factor(genotype), y=(log10(rnabl)))) + geom_boxplot() + geom_point() + theme_minimal()

#combine all pvalues into a signle tibble
getwd()
ESPRIT_pandqvalues <- left_join(tidyres_cd4_bh, tidyres_rna_bh, by = "SNPs") %>% 
  left_join(., tidyres_il6_bh, by="SNPs") %>%  left_join(., tidyres_crp_bh, by = "SNPs") %>% 
  left_join(., tidyres_ddimer_bh, by= "SNPs")

write_tsv(ESPRIT_pandqvalues, "ESPRIT_pandqvalues_feb15")



