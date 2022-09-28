# PAC2022Murray__TET2_HIV
Supporting code for the manuscript "Association between TET2 genetic variation and viral load in people living with HIV Association between TET2 genetic variation and viral load in people living with HIV"



| Trial    | Type of analysis                   | Old name                                | New name                                       |
|----------|------------------------------------|-----------------------------------------|------------------------------------------------|
| START    | Create genotype and phenotype file | START_LME_sigsnps_followup              | TET2manuscript_START_create_phenogeno          |
|          | main GLM                           | START_TET2GLM_feb15data                 | START_TET2GLM_feb15data                        |
|          | Sensitivity LME                    | START_LME_sigsnps_followup              | TET2manuscript_START_LME_sensitivity           |
|          | sensitivity black only             | START_TET2GLM_feb15data_blackonly       | TET2manuscript_START_GLM_sensitivity_blackonly |
|          | sensitivity cryptic relatedness    | START_TET2GLM_feb15data_nocryptic       |                                                |
|          | sensitivity recent infection       | START_TET2GLM_feb15data_recentinfection | TET2manuscript_START_GLM_sensitivity_recentinf |
|          | Merging with rsID and presentation |                                         |                                                |
|          | LD                                 | START_sigsnps_LD                        | TET2manuscript_START_LD                        |
| FIRST    | Create genotype and phenotype file | FIRST_create_phenogeno_feb15data        | TET2manuscript_FIRST_create_phenogeno          |
|          | main GLM                           | FIRST_TET2GLM_feb15data                 | TET2manuscript_FIRST_GLM_main                  |
|          | sensitivity black only             | FIRST_TET2GLM_feb15data_blackonly       | TET2manuscript_FIRST_GLM_sensitivity_blackonly |
|          | Merging with rsID and presentation |                                         |                                                |
|          | LD                                 | FIRST_sigsnps_LD                        |                                                |
| SMART    | Create genotype and phenotype file | SMART_create_phenogeno_feb15data        | TET2manuscript_SMART_create_phenogeno          |
|          | main GLM                           | SMART_TET2GLM_feb15                     | TET2manuscript_SMART_GLM_main                  |
| ESPRIT   | Create genotype and phenotype file |                                         |                                                |
|          | main GLM                           | ESPRIT_TET2GLM_feb15                    | TET2manuscript_ESPRIT_GLM_main                 |
| STALWART | Create genotype and phenotype file | STALWART_create_phenogeno_feb15         | TET2manuscript_STALWART_create_phenogeno       |
|          | main GLM                           | STALWART_TET2GLM_feb15                  |                                                |