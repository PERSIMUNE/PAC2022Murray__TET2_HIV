
setwd("/home/projects/cu_10103/scratch/TET2CGA/START")

library(tidyverse)
library(readxl)
library(broom)


### Load in genotype and phenotype data
START_genopheno <- read_tsv("START_genopheno_15febdata.tsv")
colnames(START_genopheno)

### Create format for multiple linear regression

## selecting and grouping SNPs
linregtet2 <-
  START_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs)
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)

## formatting data types
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)

## visually inspect result
table(linregtet2$genotype)



## check association between groups aka SNPs and variable (VL)

results_rna <-
  linregtet2 %>% do(fitSNPs = glm(log10(RNA_00) ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_rna <- tidy(results_rna, fitSNPs)

tidyres_rna_man <- tidyres_rna %>%
  filter(term == "genotype") %>%
  select("SNPs", "estimate", "std.error", "p.value") %>%
  rename(
    rna_pvalue = p.value,
    rna_estimate = estimate,
    rna_stderror = std.error
  )

tidyres_rna_bh <-
  tidyres_rna_man %>%
  ungroup(tidyres_rna_man$rna_palue) %>%
  mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method = "BH"))




## check association between groups aka SNPs and variable (CD4)

results_cd4 <-
  linregtet2 %>% do(fitSNPs = glm(CD4_00 ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_cd4 <- tidy(results_cd4, fitSNPs)

tidyres_cd4_man <- tidyres_cd4 %>%
  filter(term == "genotype") %>%
  select("SNPs", "estimate", "std.error", "p.value") %>%
  rename(
    cd4_pvalue = p.value,
    cd4_estimate = estimate,
    cd4_stderror = std.error
  )

tidyres_cd4_bh <-
  tidyres_cd4_man %>%
  ungroup(tidyres_cd4_man$cd4_palue) %>%
  mutate(cd4_qvalue = p.adjust(tidyres_cd4_man$cd4_pvalue, method = "BH"))

tidyres_cd4_bon <-
  tidyres_cd4 %>%
  filter(p.value < 1.72e-4) %>%
  filter(term == "genotype")

## check association between groups aka SNPs and variable (CD8)

results_cd8 <-
  linregtet2 %>% do(fitSNPs = glm(CD8_00 ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_cd8 <- tidy(results_cd8, fitSNPs)

tidyres_cd8_man <- tidyres_cd8 %>%
  filter(term == "genotype") %>%
  select("SNPs", "estimate", "std.error", "p.value") %>%
  rename(
    cd8_pvalue = p.value,
    cd8_estimate = estimate,
    cd8_stderror = std.error
  )

tidyres_cd8_bh <-
  tidyres_cd8_man %>%
  ungroup(tidyres_cd8_man$cd8_palue) %>%
  mutate(cd8_qvalue = p.adjust(tidyres_cd8_man$cd8_pvalue, method = "BH"))



## check association between groups aka SNPs and variable (CD4 nadir)

results_cd4nadir <-
  linregtet2 %>% do(fitSNPs = glm(cd4nadir_00 ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_cd4nadir <- tidy(results_cd4nadir, fitSNPs)

tidyres_cd4nadir_man <- tidyres_cd4nadir %>%
  filter(term == "genotype") %>%
  select("SNPs", "estimate", "std.error", "p.value") %>%
  rename(
    cd4nadir_pvalue = p.value,
    cd4nadir_estimate = estimate,
    cd4nadir_stderror = std.error
  )

tidyres_cd4nadir_bh <-
  tidyres_cd4nadir_man %>%
  ungroup(tidyres_cd4nadir_man$cd4nadir_palue) %>%
  mutate(cd4nadir_qvalue = p.adjust(tidyres_cd4nadir_man$cd4nadir_pvalue,
    method =
      "BH"
  ))


## check association between groups aka SNPs and variable (CD4/CD8ratio)

results_cd4cd8 <-
  linregtet2 %>% do(fitSNPs = glm(cd4cd8_00 ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_cd4cd8 <- tidy(results_cd4cd8, fitSNPs)

tidyres_cd4cd8_man <- tidyres_cd4cd8 %>%
  filter(term == "genotype") %>%
  select("SNPs", "estimate", "std.error", "p.value") %>%
  rename(
    cd4cd8_pvalue = p.value,
    cd4cd8_estimate = estimate,
    cd4cd8_stderror = std.error
  )

tidyres_cd4cd8_bh <-
  tidyres_cd4cd8_man %>%
  ungroup(tidyres_cd4cd8_man$cd4cd8_palue) %>%
  mutate(cd4cd8_qvalue = p.adjust(tidyres_cd4cd8_man$cd4cd8_pvalue,
    method =
      "BH"
  ))

## check association between groups aka SNPs and variable (IL6)

results_il6 <-
  linregtet2 %>% do(fitSNPs = glm(log2(IL6) ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_il6 <- tidy(results_il6, fitSNPs)

tidyres_il6_man <- tidyres_il6 %>%
  filter(term == "genotype") %>%
  select("SNPs", "estimate", "std.error", "p.value") %>%
  rename(
    il6_pvalue = p.value,
    il6_estimate = estimate,
    il6_stderror = std.error
  )

tidyres_il6_bh <-
  tidyres_il6_man %>%
  ungroup(tidyres_il6_man$il6_palue) %>%
  mutate(il6_qvalue = p.adjust(tidyres_il6_man$il6_pvalue, method = "BH"))

tidyres_il6_bon <-
  tidyres_il6 %>%
  filter(p.value < 1.72e-4) %>%
  filter(term == "genotype")

## check association between groups aka SNPs and variable (crp)

results_crp <-
  linregtet2 %>% do(fitSNPs = glm(log2(crp) ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_crp <- tidy(results_crp, fitSNPs)

tidyres_crp_man <- tidyres_crp %>%
  filter(term == "genotype") %>%
  select("SNPs", "estimate", "std.error", "p.value") %>%
  rename(
    crp_pvalue = p.value,
    crp_estimate = estimate,
    crp_stderror = std.error
  )

tidyres_crp_bh <-
  tidyres_crp_man %>%
  ungroup(tidyres_crp_man$crp_palue) %>%
  mutate(crp_qvalue = p.adjust(tidyres_crp_man$crp_pvalue, method = "BH"))


## check association between groups aka SNPs and variable (ddimer)

results_ddimer <-
  linregtet2 %>% do(fitSNPs = glm(log2(ddimer) ~ genotype + PC1 + PC2 + PC3 +
    PC4 + GENDER, data = .))
tidyres_ddimer <- tidy(results_ddimer, fitSNPs)

tidyres_ddimer_man <- tidyres_ddimer %>%
  filter(term == "genotype") %>%
  select("SNPs", "p.value") %>%
  rename(ddimer_pvalue = p.value)

tidyres_ddimer_bh <-
  tidyres_ddimer_man %>%
  ungroup(tidyres_ddimer_man$ddimer_palue) %>%
  mutate(ddimer_qvalue = p.adjust(tidyres_ddimer_man$ddimer_pvalue,
    method =
      "BH"
  ))




# combine all results into a signle tibble
getwd()
START_pandqvalues <-
  left_join(tidyres_cd4_bh, tidyres_cd4cd8_bh, by = "SNPs") %>%
  left_join(., tidyres_cd4nadir_bh,
    by =
      "SNPs"
  ) %>%
  left_join(., tidyres_cd8_bh, by = "SNPs") %>%
  left_join(., tidyres_rna_bh, by = "SNPs") %>%
  left_join(., tidyres_il6_bh, by = "SNPs") %>%
  left_join(., tidyres_crp_bh, by = "SNPs") %>%
  left_join(., tidyres_ddimer_bh, by = "SNPs")

write_tsv(START_pandqvalues, "START_pandqvalues.tsv")
