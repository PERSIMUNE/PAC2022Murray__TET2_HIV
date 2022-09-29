

library(tidyverse)
library(readxl)
library(broom)


setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
START_genopheno  <- read_tsv("START_sigsnps_cancer_genopheno.tsv")
colnames(START_genopheno)

###select black participants

START_genopheno <- START_genopheno %>% filter(START_genopheno$race == 1)

linregtet2 <-
  START_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs)
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)
table(linregtet2$genotype)



##check association between groups aka SNPs and variable (VL)


results_rna <-
  linregtet2 %>% do(fitSNPs = glm(log10(RNA_00) ~ genotype +
                                    PC1 + PC2 + PC3 +   PC4 + 
                                   GENDER.x
                                  , data = .))
tidyres_rna <- tidy(results_rna, fitSNPs)

tidyres_rna_man <- tidyres_rna %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% rename(
    rna_pvalue = p.value,
    rna_estimate = estimate,
    rna_stderror = std.error
  )

tidyres_rna_bh <-
  tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue) %>%
  mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method = "BH"))

tidyres_rna_bon <-
  tidyres_rna %>% filter(p.value < 1.72e-4) %>% filter(term == "genotype")

###save
write_tsv(tidyres_rna_bh, "START_blackonly_pandqvalues.tsv")
