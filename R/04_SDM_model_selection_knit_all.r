# 04_SDM_model_selection_knit_all

# This script knits all of the R markdown scripts in this folder, to conduct the model selection process for each taxon.

library(rmarkdown)
library(here)

# 04.1_csyif_SDM
render(input = here::here("R", "04_SDM_model_selection", "04.1_csyif_SDM", "04.1_csyif_SDM.rmd"),
       knit_root_dir = here::here())

# 04.2_cssif_SDM
render(input = here::here("R", "04_SDM_model_selection", "04.2_cssif_SDM", "04.2_cssif_SDM.rmd"),
       knit_root_dir = here::here())

# 04.3_rf_SDM
render(input = here::here("R", "04_SDM_model_selection", "04.3_rf_SDM", "04.3_rockfish_SDM.rmd"),
       knit_root_dir = here::here())

# 04.4_shrimp_larvae_SDM
render(input = here::here("R", "04_SDM_model_selection", "04.4_shrimp_larvae_SDM", "04.4_shrimp_larvae_SDM.rmd"),
       knit_root_dir = here::here())

# 04.5_cancer_crab_larvae_SDM
render(input = here::here("R", "04_SDM_model_selection", "04.5_cancer_crab_larvae_SDM", "04.5_cancer_crab_larvae_SDM.rmd"),
       knit_root_dir = here::here())

# 04.6_non_cancer_crab_larvae_SDM
render(input = here::here("R", "04_SDM_model_selection", "04.6_non_cancer_crab_larvae_SDM", "04.6_non_cancer_crab_larvae_SDM.rmd"),
       knit_root_dir = here::here())

# 04.7_hyperiid_amphipod
render(input = here::here("R", "04_SDM_model_selection", "04.7_hyperiid_amphipod_SDM", "04.7_hyperiid_amphipod_SDM.rmd"),
       knit_root_dir = here::here())

# 04.8_common_murre
render(input = here::here("R", "04_SDM_model_selection", "04.8_common_murre_SDM", "04.8_common_murre_SDM.rmd"),
       knit_root_dir = here::here())

# 04.9_sooty_shearwater
render(input = here::here("R", "04_SDM_model_selection", "04.9_sooty_shearwater_SDM", "04.9_sooty_shearwater_SDM.rmd"),
       knit_root_dir = here::here())

# 04.10_hake
render(input = here::here("R", "04_SDM_model_selection", "04.10_hake_SDM", "04.10_hake_SDM.rmd"),
       knit_root_dir = here::here())
