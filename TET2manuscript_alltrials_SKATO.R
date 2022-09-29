# Originally Written by Cameron MacPherson
# Modified by Preston Leung (Take in multiple entrez ID rather than just a single one)

################
# set work dir #
################
setwd("PATH/TO/WORKDIR") # User to set work dir

#- Load functions
source("cgahelper.R")

#- Define map function | retuns a scalar
update.cohort = function(cohort) last.cohort <<- cohort
reset.cohort = function() last.cohort <<- "NONE"

do = function(entrez, window, pidfile, confounders, cohort, outcome, ignore = 0) {
  #- Inform user
  cat(
    sprintf(
      "> Cohort=%s, Outcome=%s, Entrez=%s, Window=%s, Confounders=%s.\n",
      cohort, 
      outcome, 
      entrez, 
      window, 
      paste0("[",paste0(confounders, collapse="+"),"]")
    )
  )
  message("Step 1 - Load data")
  #- Load data
  if (cohort != last.cohort) {
    cat(sprintf("loading %s...", cohort))
    loadData(cohort=cohort)
    cat("LOADED!\n")
  }
  update.cohort(cohort)
  setDefault.TestMutator(log) # log is the log function
  
  #- Exit early if outcome not in cohort (the column title is specific to each cohort data)
  # e.g. START cohort outcome of rnabase line is named RNA_00. 
  # while in FIRST it is named rnabl
  if (!(outcome %in% colnames(DATA__$PHENOTYPES))) {
    cat("Outcome not found in cohort, returning 2.\n")
    return (2)
  }
  
  #- Specify query. pids are patient ids from the cohorts 
  pids     = read.table(pidfile)[,1]
  
  message("Step 2 - Setup test")  
  #- Setup test - here might be able to modify the entrez param for 1+ genes
  setupTest.withEntrez(
    # entrez = entrez, 
    entrez = unlist(stringr::str_split(entrez,"_")),
    controls = confounders, 
    outcome = outcome, 
    pad = window, 
    pids = pids
  )
  
  TEST__$y[is.infinite(TEST__$y)] <<- NA
  
  message("Step 3 - Run SKATO test")
  #- Run test
  runSKAT.SKATO()
  #- Return p.value
  result = getResultsSKAT.SKATO()$p.value
  #- Inform user
  cat(sprintf("p-value = %s", result))
  cat("=========\n\n")
  result
}


####################################################
#- Read in parameters (Change parameter file here) #
####################################################
parameters = read.delim("/PATH/TO/Analysis.TET2.parameters.csv",sep=",",as.is=TRUE) # User to provide parameter file

#- Tidy magic to collapse confounders to a vector
parameters$confounders = parameters %>%
  dplyr::select(dplyr::starts_with("confounder.")) %>%
  purrr::pmap(function (...) {list(...)}) %>%
  purrr::map(function (listOfConfounders) as.character(listOfConfounders))

#- Remove all confounder.# columns
parameters %<>% dplyr::select(-dplyr::starts_with("confounder."))

#- Map parameters
reset.cohort()
results = purrr::pmap_dbl(parameters, do)


#- The output in 'results' is just a list of p-values without context
#- Here we make nice formatted table to match p-values to the test conducted
#- Note for Daniel: some manual formatting required afterwards to make it 
#- look like the table for the paper
result_formatted = tibble(
  entrez = parameters$entrez, 
  cohort = parameters$cohort, 
  outcomes = parameters$outcome, 
  pvalue = results
)

################################
# Change output file name here #
################################
data.table::fwrite(
  result_formatted,
  file = "Entrez_association_Outcome.tsv", 
  sep = "\t"
)


