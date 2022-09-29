library(tidyverse)
library(broom)

###To remove pariticpants with cryptic relatedness and outlying heteroxygosity, load a list of PIDs that were used in a sensitivity analysis for Ekenberg et al JID
###https://doi.org/10.1093/infdis/jiz294

###this analysis removed the 103 patients that exhibited cryptic relatedness and outlying heterozygosity + 3 additional patients that did not have RNA at baseline

setwd("/home/projects/cu_10103/scratch/results_INSIGHT_new_20180912/START/phenotypic_data")

restrictive <- read_tsv("hla_START_retrictive_2440_Feb15.tsv") %>% select(1)

vec <- rep(1, length = 2440)

restrictive <- cbind(restrictive, vec)
restrictive <- restrictive %>% rename("restrictive" = "vec")


###tet2 patient list - all patients from my TET2 proejct

setwd("/home/projects/cu_10103/scratch/TET2CGA/START")

tet <- read_tsv("START_genopheno_15febdata.tsv")%>% select(1, "RNA_00")
print(tet$RNA_00
      )
###remove the three patients without RNA (as these were also removed christinas paper)
tet$ID <- as.character(tet$ID)
restrictive$ID <- as.character(restrictive$ID)
tet <- tet %>% filter(RNA_00 > 0)

vec2 <- rep(1, length = 2543)
tet <- cbind(tet, vec2)

cryptic <- anti_join(tet, restrictive, by = "ID") %>% rename("cryptic" = "vec2") %>% select(-RNA_00)



setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
START_genopheno  <- read_tsv("START_sigsnps_cancer_genopheno.tsv")
colnames(START_genopheno)

START_genopheno$ID <- as.character(START_genopheno$ID)
START_genopheno <- left_join(START_genopheno, cryptic, by = "ID") 

START_genopheno$cryptic <- START_genopheno$cryptic %>% replace_na(0)

print(START_genopheno$cryptic)

START_genopheno <- START_genopheno %>% filter(cryptic == 0)

###Visually check numbers match up based


#set up for GLM

linregtet2 <- START_genopheno %>% gather("SNPs", "genotype", starts_with("AX-"))
linregtet2 <- linregtet2 %>% group_by(SNPs)
# linregtet2 <- linregtet2 %>% mutate(genotype = factor(genotype)) %>% group_by(SNPs) 
# levels(linregtet2$genotype) <- list("0" = 0, "1" = 1, "2" = 2, "NA" = NA)
levels(linregtet2$genotype)
linregtet2$genotype <- as.integer(linregtet2$genotype)
table(linregtet2$genotype)

##check association between groups aka SNPs and variable (VL) - main model

results_rna <- linregtet2 %>% do(fitSNPs = glm(log10(RNA_00) ~ genotype + 
                                                 PC1 + PC2 + PC3 +PC4 + 
                                                 GENDER.x, data = .))
tidyres_rna <- tidy(results_rna, fitSNPs)

tidyres_rna_man <- tidyres_rna %>%
  filter(term == "genotype") %>% select("SNPs", "estimate", "std.error", "p.value") %>% 
  rename(rna_pvalue = p.value, rna_estimate = estimate, rna_stderror = std.error)


tidyres_rna_bh <- tidyres_rna_man %>% ungroup(tidyres_rna_man$rna_palue) %>% 
  mutate(rna_qvalue = p.adjust(tidyres_rna_man$rna_pvalue, method="BH"))


write_tsv(tidyres_rna_bh, "START_nocryptic_pandqvalues.tsv")


