setwd("/home/projects/cu_10103/scratch/TET2CGA/SMART")
library(tidyverse)
library(readxl)
library(broom)

SMART_genopheno <- read_tsv("SMART_genopheno_feb15.tsv")

linregtet2 <- SMART_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs)
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)
table(linregtet2$genotype)

## check association between groups aka SNPs and variable (VL)
### Remove as SMART patients are not ART naive
# results_rna <- linregtet2 %>% do(fitSNPs = glm(log10(rnabl) ~ genotype + PC1 + PC2 + PC3 +PC4 + GENDER, data = .))
# tidyres_rna <- tidy(results_rna, fitSNPs)

# tidyres_rna_man <- tidyres_rna %>%
filter(term == "genotype") %>%
  select("SNPs", "p.value") %>%
  rename(rna_pvalue = p.value)

# tidyres_rna_bh <- tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue) %>%
# mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method="BH"))


## check association between groups aka SNPs and variable (IL6)

results_il6 <- linregtet2 %>% do(fitSNPs = glm(log2(il6bl) ~ genotype + PC1 + PC2 + PC3 + PC4 + GENDER, data = .))
tidyres_il6 <- tidy(results_il6, fitSNPs)

tidyres_il6_man <- tidyres_il6 %>%
  filter(term == "genotype") %>%
  select("SNPs", "p.value") %>%
  rename(il6_pvalue = p.value)

tidyres_il6_bh <- tidyres_il6_man %>%
  ungroup(tidyres_il6_man$il6_palue) %>%
  mutate(il6_qvalue = p.adjust(tidyres_il6_man$il6_pvalue, method = "BH"))


## check association between groups aka SNPs and variable (CD4)

results_cd4 <- linregtet2 %>% do(fitSNPs = glm(cd4bl ~ genotype + PC1 + PC2 + PC3 + PC4 + GENDER, data = .))
tidyres_cd4 <- tidy(results_cd4, fitSNPs)

tidyres_cd4_man <- tidyres_cd4 %>%
  filter(term == "genotype") %>%
  select("SNPs", "p.value") %>%
  rename(cd4_pvalue = p.value)

tidyres_cd4_bh <- tidyres_cd4_man %>%
  ungroup(tidyres_cd4_man$cd4_palue) %>%
  mutate(cd4_qvalue = p.adjust(tidyres_cd4_man$cd4_pvalue, method = "BH"))


## check association between groups aka SNPs and variable (crp)

results_crp <- linregtet2 %>% do(fitSNPs = glm(log2(crpbl) ~ genotype + PC1 + PC2 + PC3 + PC4 + GENDER, data = .))
tidyres_crp <- tidy(results_crp, fitSNPs)

tidyres_crp_man <- tidyres_crp %>%
  filter(term == "genotype") %>%
  select("SNPs", "p.value") %>%
  rename(crp_pvalue = p.value)

tidyres_crp_bh <- tidyres_crp_man %>%
  ungroup(tidyres_crp_man$crp_palue) %>%
  mutate(crp_qvalue = p.adjust(tidyres_crp_man$crp_pvalue, method = "BH"))


## check association between groups aka SNPs and variable (ddimer)

results_ddimer <- linregtet2 %>% do(fitSNPs = glm(log2(ddimerbl) ~ genotype + PC1 + PC2 + PC3 + PC4 + GENDER, data = .))
tidyres_ddimer <- tidy(results_ddimer, fitSNPs)

tidyres_ddimer_man <- tidyres_ddimer %>%
  filter(term == "genotype") %>%
  select("SNPs", "p.value") %>%
  rename(ddimer_pvalue = p.value)

tidyres_ddimer_bh <- tidyres_ddimer_man %>%
  ungroup(tidyres_ddimer_man$ddimer_palue) %>%
  mutate(ddimer_qvalue = p.adjust(tidyres_ddimer_man$ddimer_pvalue, method = "BH"))


## check association between groups aka SNPs and variable (cd8) - note, lots of missing data for cd8 and so data not presented

results_cd8 <- linregtet2 %>% do(fitSNPs = glm(cd8 ~ genotype + PC1 + PC2 + PC3 + PC4 + GENDER, data = .))
tidyres_cd8 <- tidy(results_cd8, fitSNPs)

tidyres_cd8_man <- tidyres_cd8 %>%
  filter(term == "genotype") %>%
  select("SNPs", "p.value") %>%
  rename(cd8_pvalue = p.value)

tidyres_cd8_bh <- tidyres_cd8_man %>%
  ungroup(tidyres_cd8_man$cd8_palue) %>%
  mutate(cd8_qvalue = p.adjust(tidyres_cd8_man$cd8_pvalue, method = "BH"))


# plot

linregtet2 %>%
  filter(SNPs == "AX-14431780") %>%
  ggplot(aes(x = factor(genotype), y = (log10(rnabl)))) +
  geom_boxplot() +
  geom_point() +
  theme_minimal()

# combine all pvalues into a signle tibble
getwd()
SMART_pandqvalues <- left_join(tidyres_cd4_bh, tidyres_cd8_bh, by = "SNPs") %>%
  left_join(., tidyres_rna_bh, by = "SNPs") %>%
  left_join(., tidyres_il6_bh, by = "SNPs") %>%
  left_join(., tidyres_crp_bh, by = "SNPs") %>%
  left_join(., tidyres_ddimer_bh, by = "SNPs")

write_tsv(SMART_pandqvalues, "SMART_pandqvalues_feb15")
