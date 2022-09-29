library(tidyverse)
library(lubridate)
library(nlme)

### import main START phenotype and include all RNA values, dates and censor visits depending on whether particpiants have started therapy at particular visits
setwd("/home/projects/cu_10103/data/INSIGHT_genetics/START/START_27102017")

# Import START phenotype data including
STARTgenomics <- read_tsv("startgenomics.tsv")

### only select relevant data
START_basic <-
  select(
    STARTgenomics,
    "ID",
    "GENDER",
    "age_00",
    starts_with("lrna_"),
    starts_with("RNADT_"),
    "STARTDT",
    "race",
    "rtreat",
    "CD4_00",
    "cd4cd8_00",
    "cgroup",
  )
### load START sigSNPs
setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
start_geno <- read_tsv("START_sigsnps_cancer_genopheno.tsv") %>% select("ID", starts_with("AX"), starts_with("PC"))


### merge

START_LME_sigsnps <- left_join(start_geno, START_basic, by = "ID")

# Turn from wide to long DF - based on visits - this is to set up the DF in a way so I can censor  visits where the participant has already started therapy

STARTlong <- START_LME_sigsnps %>% gather(visit, value, starts_with("lrna"))
STARTlong$STARTDT <- as_date(mdy((STARTlong$STARTDT)))


### create individual dataframes for each visit to create the days between ART and RNA variable for each visit

####### baseline############

STARTlong_00 <- STARTlong %>% filter(visit == "lrna_00")

### create variable indicating how long a visit is
STARTlong_00$RNADT_00 <- as_date(mdy(STARTlong_00$RNADT_00))


STARTlong_00 <- STARTlong_00 %>% mutate("daysbetweenARTandRNA" = (interval(start = STARTlong_00$RNADT_00, end = STARTlong_00$STARTDT) / duration(num = 1, units = "days")))

### DO not filter visits for Baseline, as it is just removing all patients with NA start date

# STARTlong_00 <- STARTlong_00 %>% filter(daysbetweenARTandRNA > 0)

############ month 1###################



STARTlong_01 <- STARTlong %>% filter(visit == "lrna_01")

### create variable indicating how long a visit is
STARTlong_01$RNADT_01 <- as_date(mdy(STARTlong_01$RNADT_01))


STARTlong_01 <- STARTlong_01 %>% mutate("daysbetweenARTandRNA" = (interval(start = STARTlong_01$RNADT_01, end = STARTlong_01$STARTDT) / duration(num = 1, units = "days")))

### filter for visits that have not started yet

STARTlong_01 <- STARTlong_01 %>% filter(daysbetweenARTandRNA > 0)


##### Month 4############


STARTlong_04 <- STARTlong %>% filter(visit == "lrna_04")

### create variable indicating how long a visit is
STARTlong_04$RNADT_04 <- as_date(mdy(STARTlong_04$RNADT_04))


STARTlong_04 <- STARTlong_04 %>% mutate("daysbetweenARTandRNA" = (interval(start = STARTlong_04$RNADT_04, end = STARTlong_04$STARTDT) / duration(num = 1, units = "days")))

### filter for visits that have not started yet

STARTlong_04 <- STARTlong_04 %>% filter(daysbetweenARTandRNA > 0)


##### Month 8############


STARTlong_08 <- STARTlong %>% filter(visit == "lrna_08")

### create variable indicating how long a visit is
STARTlong_08$RNADT_08 <- as_date(mdy(STARTlong_08$RNADT_08))


STARTlong_08 <- STARTlong_08 %>% mutate("daysbetweenARTandRNA" = (interval(start = STARTlong_08$RNADT_08, end = STARTlong_08$STARTDT) / duration(num = 1, units = "days")))

### filter for visits that have not started yet

STARTlong_08 <- STARTlong_08 %>% filter(daysbetweenARTandRNA > 0)

##### Month 12############


STARTlong_12 <- STARTlong %>% filter(visit == "lrna_12")

### create variable indicating how long a visit is
STARTlong_12$RNADT_12 <- as_date(mdy(STARTlong_12$RNADT_12))


STARTlong_12 <- STARTlong_12 %>% mutate("daysbetweenARTandRNA" = (interval(start = STARTlong_12$RNADT_12, end = STARTlong_12$STARTDT) / duration(num = 1, units = "days")))

### filter for visits that have not started yet

STARTlong_12 <- STARTlong_12 %>% filter(daysbetweenARTandRNA > 0)

### remove dates to bind it together (something about the date formating in the original dataframes make)
STARTlong_00 <- STARTlong_00 %>% select(-starts_with("RNADT_"))
STARTlong_01 <- STARTlong_01 %>% select(-starts_with("RNADT_"))
STARTlong_04 <- STARTlong_04 %>% select(-starts_with("RNADT_"))
STARTlong_08 <- STARTlong_08 %>% select(-starts_with("RNADT_"))
STARTlong_12 <- STARTlong_12 %>% select(-starts_with("RNADT_"))

####### Rbind it all together

START_long_all <- bind_rows(STARTlong_00, STARTlong_01, STARTlong_04, STARTlong_08, STARTlong_12)


######## Run LME#############
### convert to dataframe
START_long_all <- as.data.frame(START_long_all)

## remove the -'s in the data frame as these arent compatable with the lme model

names(START_long_all) <- gsub(x = names(START_long_all), pattern = "\\-", replacement = ".")

### run lme on a single SNP to get structure and check its working. Summarise results in a new dataframe

out1 <- data.frame(NULL) # create empty dataframe to store results
m <- lme(value ~ factor(visit) + PC1 + PC2 + PC3 + PC4 + GENDER + AX.14431677, random = ~ 1 | ID, na.action = na.omit, data = START_long_all)

summary(m)
sum <- summary(m)
tidy <- broom::tidy(m)
sumTtable <- sum$tTable


out1[1, 1] <- "AX.14431677" # print variable name
out1[1, 2] <- sum$tTable[11, 1] # intercept
out1[1, 3] <- sum$tTable[11, 2] # standard error
out1[1, 4] <- sum$tTable[11, 5] # P value

###Run on all SNPs using a loop loop https://stackoverflow.com/questions/41389799/loops-in-r-linear-regression
out <- data.frame(SNP = character(33), Estimate = numeric(33), Std.error = numeric(33), p.value = numeric(33)) # create object to store resutls

snps <- names(START_long_all[c(2:34)])
print(snps)

for (i in 1:length(snps)) {
  START_long_all$SNP <- as.numeric(unlist(START_long_all[, snps[i]]))

  m <- lme(value ~ factor(visit) + PC1 + PC2 + PC3 + PC4 + GENDER + SNP, random = ~ 1 | ID, na.action = na.exclude, data = START_long_all)
  m_sum <- summary(m)
  out[i, 1] <- snps[i]
  print(snps[i]) # print variable name
  out[i, 2] <- m_sum$tTable[11, 1] # intercept
  out[i, 3] <- m_sum$tTable[11, 2] # standard error
  out[i, 4] <- m_sum$tTable[11, 5] # P value
}

outround <- round(out[2:4], 5)
LME_summary <- cbind(snps, outround)

LME_summary$snps <- str_replace(LME_summary$snps, "AX.", "AX-")

setwd("/home/projects/cu_10103/scratch/TET2CGA/START")
write_tsv(LME_summary, "START_followup_sigsnips_LME.tsv")

