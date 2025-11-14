## 11_Fig3_covariate_time_series

# Description: This script generates Figure 3, which illustrates the 
# time series of the different SDM-derived marine survival covariates.


## Load libraries
library(tidyverse)
library(readxl)
library(here)
library(viridis)
library(broom)
library(ggpubr)
library(sf)
library(lubridate)

# source the prep scripts for each of the surveys
# source(here::here("R", "CCES_make_mesh.R"))
source(here::here("R", "PRS_PWCC_make_mesh.R"))
source(here::here("R", "hake_survey_make_mesh.R"))
source(here::here("R", "JSOES_seabirds_make_mesh.R"))
source(here::here("R", "JSOES_make_mesh.R"))


# load the SDM outputs + estimated uncertainty from the stage 1 models
# 05.1_prey_SDM
load(here::here("R", "05_stage1_SDM", "05.1_prey_field_SDM", "prey_field_SDM_output.rda"))
prey_field_SDM_Obj <- prey_field_SDM_output$prey_field_SDM_Obj
prey_field_SDM_Opt <- prey_field_SDM_output$prey_field_SDM_Opt
prey_field_SDM_report <- prey_field_SDM_output$prey_field_SDM_report

load(here::here("R", "05_stage1_SDM", "05.1_prey_field_SDM", "estimated_SE_prey_model.rda"))
SE_pianka_o_csyif_prey_field_t_prey_model <- SE_prey_model$SE_pianka_o_csyif_prey_field_t_prey_model
SE_pianka_o_cssif_prey_field_t_prey_model <- SE_prey_model$SE_pianka_o_cssif_prey_field_t_prey_model

# 05.2_seabird_SDM
load(here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "seabird_SDM_output.rda"))
seabird_SDM_Obj <- seabird_SDM_output$seabird_SDM_Obj
seabird_SDM_Opt <- seabird_SDM_output$seabird_SDM_Opt
seabird_SDM_report <- seabird_SDM_output$seabird_SDM_report

load(here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "estimated_SE_seabird_model.rda"))
SE_pianka_o_cssif_sosh_t_seabird_model <- SE_seabird_model$SE_pianka_o_cssif_sosh_t_seabird_model
SE_pianka_o_cssif_comu_t_seabird_model <- SE_seabird_model$SE_pianka_o_cssif_comu_t_seabird_model
SE_pianka_o_csyif_sosh_t_seabird_model <- SE_seabird_model$SE_pianka_o_csyif_sosh_t_seabird_model
SE_pianka_o_csyif_comu_t_seabird_model <- SE_seabird_model$SE_pianka_o_csyif_comu_t_seabird_model

# 05.3_hake_SDM
# load(here::here("R", "05_stage1_SDM", "05.3_hake_SDM", "hake_SDM_output.rda"))
# hake_SDM_Obj <- hake_SDM_output$hake_SDM_Obj
# hake_SDM_Opt <- hake_SDM_output$hake_SDM_Opt
# hake_SDM_report <- hake_SDM_output$hake_SDM_report
# 
# load(here::here("R", "05_stage1_SDM", "05.3_hake_SDM", "estimated_SE_hake_model.rda"))
# SE_pianka_o_csyif_hake_t_hake_model <- SE_hake_model$SE_pianka_o_csyif_hake_t_hake_model
# SE_pianka_o_cssif_hake_t_hake_model <- SE_hake_model$SE_pianka_o_cssif_hake_t_hake_model



#### identify where in our covariance matrix our ADREPORTed variables are ####

### Prey field model
common_years_prey <- intersect(
  unique(jsoes_bongo_cancer_crab_larvae$year), # jsoes bongo
  unique(rf$year) # PRS/PWCC
)

# first extract the summary of the SD
prey_field_SDM_SD_summary <- summary(prey_field_SDM_Opt$SD)

# use this object to figure out what the standard errors are from our variables of interest, and then match them to the cov object
# csyif index of abundance
csyif_index_of_abundance_indices <- which(grepl("csyif_index_of_abundance", rownames(prey_field_SDM_SD_summary)))
csyif_index_of_abundance_cov_indices <- which(sqrt(diag(prey_field_SDM_Opt$SD$cov)) %in% prey_field_SDM_SD_summary[csyif_index_of_abundance_indices, "Std. Error"])

# cssif index of abundance
cssif_index_of_abundance_indices <- which(grepl("cssif_index_of_abundance", rownames(prey_field_SDM_SD_summary)))
cssif_index_of_abundance_cov_indices <- which(sqrt(diag(prey_field_SDM_Opt$SD$cov)) %in% prey_field_SDM_SD_summary[cssif_index_of_abundance_indices, "Std. Error"])

# prey field
prey_field_index_of_abundance_indices <- which(grepl("prey_field_index_of_abundance", rownames(prey_field_SDM_SD_summary)))
prey_field_index_of_abundance_cov_indices <- which(sqrt(diag(prey_field_SDM_Opt$SD$cov)) %in% prey_field_SDM_SD_summary[prey_field_index_of_abundance_indices, "Std. Error"])

# extract the covariance matrices for each using the indices
# csyif
csyif_index_of_abundance_cov_matrix <- as.matrix(prey_field_SDM_Opt$SD$cov[csyif_index_of_abundance_cov_indices, csyif_index_of_abundance_cov_indices])
colnames(csyif_index_of_abundance_cov_matrix) <- sort(unique(csyif$year))
rownames(csyif_index_of_abundance_cov_matrix) <- sort(unique(csyif$year))

# cssif
cssif_index_of_abundance_cov_matrix <- as.matrix(prey_field_SDM_Opt$SD$cov[cssif_index_of_abundance_cov_indices, cssif_index_of_abundance_cov_indices])
colnames(cssif_index_of_abundance_cov_matrix) <- sort(unique(cssif$year))
rownames(cssif_index_of_abundance_cov_matrix) <- sort(unique(cssif$year))

# prey field
prey_field_index_of_abundance_cov_matrix <- as.matrix(prey_field_SDM_Opt$SD$cov[prey_field_index_of_abundance_cov_indices, prey_field_index_of_abundance_cov_indices])
colnames(prey_field_index_of_abundance_cov_matrix) <- common_years_prey
rownames(prey_field_index_of_abundance_cov_matrix) <- common_years_prey

### Seabird model

# first extract the summary of the SD
seabird_SDM_SD_summary <- summary(seabird_SDM_Opt$SD)

# use this object to figure out what the standard errors are from our variables of interest, and then match them to the cov object
# comu index of abundance
comu_index_of_abundance_indices <- which(grepl("comu_index_of_abundance", rownames(seabird_SDM_SD_summary)))
comu_index_of_abundance_cov_indices <- which(sqrt(diag(seabird_SDM_Opt$SD$cov)) %in% seabird_SDM_SD_summary[comu_index_of_abundance_indices, "Std. Error"])

# sosh index of abundance
sosh_index_of_abundance_indices <- which(grepl("sosh_index_of_abundance", rownames(seabird_SDM_SD_summary)))
sosh_index_of_abundance_cov_indices <- which(sqrt(diag(seabird_SDM_Opt$SD$cov)) %in% seabird_SDM_SD_summary[sosh_index_of_abundance_indices, "Std. Error"])

# extract the covariance matrices for each using the indices
# comu
comu_index_of_abundance_cov_matrix <- as.matrix(seabird_SDM_Opt$SD$cov[comu_index_of_abundance_cov_indices, comu_index_of_abundance_cov_indices])
colnames(comu_index_of_abundance_cov_matrix) <- sort(unique(comu$year))
rownames(comu_index_of_abundance_cov_matrix) <- sort(unique(comu$year))

# sosh
sosh_index_of_abundance_cov_matrix <- as.matrix(seabird_SDM_Opt$SD$cov[sosh_index_of_abundance_cov_indices, sosh_index_of_abundance_cov_indices])
colnames(sosh_index_of_abundance_cov_matrix) <- min(sosh$year):max(sosh$year)
rownames(sosh_index_of_abundance_cov_matrix) <- min(sosh$year):max(sosh$year)


### Hake model

# # first extract the summary of the SD
# hake_SDM_SD_summary <- summary(hake_SDM_Opt$SD)
# 
# # use this object to figure out what the standard errors are from our variables of interest, and then match them to the cov object
# hake_index_of_abundance_indices <- which(grepl("hake_index_of_abundance", rownames(hake_SDM_SD_summary)))
# hake_index_of_abundance_cov_indices <- which(sqrt(diag(hake_SDM_Opt$SD$cov)) %in% hake_SDM_SD_summary[hake_index_of_abundance_indices, "Std. Error"])

# extract the covariance matrices for hake
# hake_index_of_abundance_cov_matrix <- as.matrix(hake_SDM_Opt$SD$cov[hake_index_of_abundance_cov_indices, hake_index_of_abundance_cov_indices])
# colnames(hake_index_of_abundance_cov_matrix) <- min(hake$year):max(hake$year)
# rownames(hake_index_of_abundance_cov_matrix) <- min(hake$year):max(hake$year)

#### Collate marine survival covariates and their uncertainties ####

### Indices of abundance
prey_field_index_df <- data.frame(year = common_years_prey,
                            prey_field_index = prey_field_SDM_report$prey_field_index_of_abundance/max(prey_field_SDM_report$prey_field_index_of_abundance),
                            prey_field_index_SE = sqrt(diag(prey_field_index_of_abundance_cov_matrix))/(max(prey_field_SDM_report$prey_field_index_of_abundance)))

csyif_index_df <- data.frame(year = sort(unique(csyif$year)),
                             csyif_index = prey_field_SDM_report$csyif_index_of_abundance/max(prey_field_SDM_report$csyif_index_of_abundance),
                             csyif_index_SE = sqrt(diag(csyif_index_of_abundance_cov_matrix))/max(prey_field_SDM_report$csyif_index_of_abundance))

cssif_index_df <- data.frame(year = sort(unique(cssif$year)),
                             cssif_index = prey_field_SDM_report$cssif_index_of_abundance/max(prey_field_SDM_report$cssif_index_of_abundance),
                             cssif_index_SE = sqrt(diag(cssif_index_of_abundance_cov_matrix))/max(prey_field_SDM_report$cssif_index_of_abundance))

comu_index_df <- data.frame(year = sort(unique(comu$year)),
                            comu_index = seabird_SDM_report$comu_index_of_abundance/max(seabird_SDM_report$comu_index_of_abundance),
                            comu_index_SE = sqrt(diag(comu_index_of_abundance_cov_matrix))/max(seabird_SDM_report$comu_index_of_abundance))

sosh_index_df <- data.frame(year = min(sosh$year):max(sosh$year),
                            sosh_index = seabird_SDM_report$sosh_index_of_abundance/max(seabird_SDM_report$sosh_index_of_abundance),
                            sosh_index_SE = sqrt(diag(sosh_index_of_abundance_cov_matrix))/max(seabird_SDM_report$sosh_index_of_abundance))

# combine these together

csyif_index_df %>% 
  left_join(cssif_index_df, by = "year") %>% 
  left_join(prey_field_index_df, by = "year") %>% 
  left_join(comu_index_df, by = "year") %>% 
  left_join(sosh_index_df, by = "year") -> indices_df

indices_df %>% 
  pivot_longer(., cols = -c("year")) %>% 
  mutate(variable = ifelse(grepl("SE", name), "SE", "estimate")) %>% 
  mutate(index = gsub("_SE", "", name)) -> indices_long

indices_long %>% 
  dplyr::select(-name) %>% 
  pivot_wider(names_from = variable, values_from = value) -> indices_long


### Overlap metrics

prey_field_csyif_overlap_df <- data.frame(year = common_years_prey,
                                          pianka_o_csyif_prey_field = prey_field_SDM_report$pianka_o_csyif_prey_field_t,
                                          pianka_o_csyif_prey_field_SE = SE_pianka_o_csyif_prey_field_t_prey_model)

comu_csyif_overlap_df <- data.frame(year = sort(unique(comu$year)),
                                          pianka_o_csyif_comu = seabird_SDM_report$pianka_o_csyif_comu_t,
                                          pianka_o_csyif_comu_SE = SE_pianka_o_csyif_comu_t_seabird_model)

sosh_csyif_overlap_df <- data.frame(year = sort(unique(sosh$year)),
                                    pianka_o_csyif_sosh = seabird_SDM_report$pianka_o_csyif_sosh_t,
                                    pianka_o_csyif_sosh_SE = SE_pianka_o_csyif_sosh_t_seabird_model)

prey_field_cssif_overlap_df <- data.frame(year = common_years_prey,
                                          pianka_o_cssif_prey_field = prey_field_SDM_report$pianka_o_cssif_prey_field_t,
                                          pianka_o_cssif_prey_field_SE = SE_pianka_o_cssif_prey_field_t_prey_model)

comu_cssif_overlap_df <- data.frame(year = sort(unique(comu$year)),
                                    pianka_o_cssif_comu = seabird_SDM_report$pianka_o_cssif_comu_t,
                                    pianka_o_cssif_comu_SE = SE_pianka_o_cssif_comu_t_seabird_model)

sosh_cssif_overlap_df <- data.frame(year = sort(unique(sosh$year)),
                                    pianka_o_cssif_sosh = seabird_SDM_report$pianka_o_cssif_sosh_t,
                                    pianka_o_cssif_sosh_SE = SE_pianka_o_cssif_sosh_t_seabird_model)

# combine these together

prey_field_csyif_overlap_df %>% 
  left_join(comu_csyif_overlap_df, by = "year") %>% 
  left_join(sosh_csyif_overlap_df, by = "year")  -> csyif_overlap_indices

prey_field_cssif_overlap_df %>% 
  left_join(comu_cssif_overlap_df, by = "year") %>% 
  left_join(sosh_cssif_overlap_df, by = "year")  -> cssif_overlap_indices

#### Generate figures ####

ggplot(indices_long, aes(x = year, y = estimate, 
                         ymin = estimate - 1.96*SE,
                         ymax = estimate + 1.96*SE)) +
  geom_point() +
  geom_errorbar() +
  geom_line() +
  facet_wrap(~index, ncol = 1)


