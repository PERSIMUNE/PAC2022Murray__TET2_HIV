

library(tidyverse)
library(readxl)
library(broom)


setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
START_genopheno  <- read_tsv("START_sigsnps_cancer_genopheno.tsv")
colnames(START_genopheno)


### import data on recent infection

setwd("/home/projects/cu_10103/data/INSIGHT_genetics/START/START_24042018")

recent <- read_tsv("recent.tsv")

###Join with START
START_genopheno <- left_join(START_genopheno,recent, by = "ID")



###remove participants that were recently infected


START_genopheno <- START_genopheno %>% filter(START_genopheno$recentinfect == 0)



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
    rna_pvalue_recent = p.value,
    rna_estimate_recent = estimate,
    rna_stderror_recent = std.error
  )

tidyres_rna_bh <-
  tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue_recent) %>%
  mutate(rna_qvalue_recent = p.adjust(tidyres_rna_man$rna_pvalue_recent, method = "BH"))



###save
setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
write_tsv(tidyres_rna_bh, "START_recentremoved_pandqvalues.tsv")
