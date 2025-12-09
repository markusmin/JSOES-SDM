## 13_EiV_fit_fig

# Description: This script generates figures showing the “fit” of the observed value 
# for covariates based on a latent covariate being estimated


## Load libraries
library(tidyverse)
library(readxl)
library(here)
library(viridis)
library(broom)
library(ggpubr)
library(sf)
library(lubridate)

# source the function scripts
source("R/functions/cAIC.R")
source("R/functions/rmvnorm_prec.R")
source("R/functions/sample_var.R")
source("R/functions/plot_anisotropy_MM.R")
source("R/functions/make_anisotropy_spde.R")
source("R/functions/H_matrix_prior_predictive_check_for_ggplot.R")
source("R/functions/plot_distribution.R")
source("R/functions/plot_distribution_PRS_PWCC.R")
source("R/functions/plot_distribution_jsoes_bongo.R")
# source("R/functions/plot_distribution_hake_survey.R")
source("R/functions/generate_prediction_maps.R")
source("R/functions/generate_prediction_maps_PRS_PWCC.R")
source("R/functions/generate_prediction_maps_jsoes_bongo.R")
source("R/functions/generate_prediction_maps_hake_survey.R")
source("R/functions/generate_prediction_maps_cces.R")
source("R/functions/stage2_helper_functions.R")


#### Load the models for Upper Columbia Summer/Fall and Snake River Fall fish run separately ####

### 06.2 seabird model

## SRF Models

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.2_seabird_SAR", "SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output_zscored.rda"))

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Obj <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Obj
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.2_seabird_SAR", "SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_output_zscored.rda"))

SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_Obj <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_Obj
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_Opt = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_Opt
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report

## UCSF Models

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.2_seabird_SAR", "UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_output_zscored.rda"))

UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_Obj <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_output$UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_Obj
UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_Opt <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_output$UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_Opt
UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_output$UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.2_seabird_SAR", "UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_output_zscored.rda"))

UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_Obj <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_output$UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_Obj
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_Opt = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_output$UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_Opt
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_output$UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report

### 06.3 hake model

## SRF Models

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.3_hake_SAR", "SRF_06_3_hake_prey_csyif_only_SAR_SEcov_output_zscored.rda"))

SRF_06_3_hake_prey_csyif_only_SAR_SEcov_Obj <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_output$SRF_06_3_hake_prey_csyif_only_SAR_SEcov_Obj
SRF_06_3_hake_prey_csyif_only_SAR_SEcov_Opt <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_output$SRF_06_3_hake_prey_csyif_only_SAR_SEcov_Opt
SRF_06_3_hake_prey_csyif_only_SAR_SEcov_report <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_output$SRF_06_3_hake_prey_csyif_only_SAR_SEcov_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.3_hake_SAR", "SRF_06_3_hake_prey_cssif_only_SAR_SEcov_output_zscored.rda"))

SRF_06_3_hake_prey_cssif_only_SAR_SEcov_Obj <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_output$SRF_06_3_hake_prey_cssif_only_SAR_SEcov_Obj
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_Opt = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_output$SRF_06_3_hake_prey_cssif_only_SAR_SEcov_Opt
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_report = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_output$SRF_06_3_hake_prey_cssif_only_SAR_SEcov_report

## UCSF Models

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.3_hake_SAR", "UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_output_zscored.rda"))

UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_Obj <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_output$UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_Obj
UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_Opt <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_output$UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_Opt
UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_report <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_output$UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.3_hake_SAR", "UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_output_zscored.rda"))

UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_Obj <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_output$UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_Obj
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_Opt = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_output$UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_Opt
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_report = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_output$UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_report


common_years_seabird_prey <- c(2003, 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2015, 2016, 2017, 2018, 2019)
jsoes_years <- c(1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 
                 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 
                 2020, 2021)

#### define a function to compare data + estimated state of a covariate ####
reformat_latent_state_df <- function(SAR_summary, SDM_summary,
                                     SAR_years, SDM_years,
                                     SAR_parameter_name,
                                     SDM_parameter_name){
  
  subset(SDM_summary, grepl(SDM_parameter_name, parameter)) -> SDM_param_summary
  subset(SAR_summary, grepl(SAR_parameter_name, parameter)) -> SAR_param_summary
  
  SDM_param_summary$year <- SDM_years
  SAR_param_summary$year <- SAR_years
  
  SDM_param_summary %>% 
    dplyr::select("Estimate", "Std. Error", "year") %>% 
    dplyr::rename(SDM_estimate = "Estimate",
                  SDM_SE = "Std. Error") -> SDM_param_summary
  
  SAR_param_summary %>% 
    dplyr::select("Estimate", "Std. Error", "year") %>% 
    dplyr::rename(latent_estimate = "Estimate",
                  latent_SE = "Std. Error") -> SAR_param_summary
  
  SDM_param_summary %>% 
    left_join(SAR_param_summary, by = "year") %>% 
    filter(!(is.na(latent_estimate))) -> eiv_comp_summary
  
  # z-score SDM estimates for comparison
  
  eiv_comp_summary$SDM_estimate_scaled = (eiv_comp_summary$SDM_estimate - mean(eiv_comp_summary$SDM_estimate))/sd(eiv_comp_summary$SDM_estimate)
  eiv_comp_summary$SDM_SE_scaled = eiv_comp_summary$SDM_SE  * (1/sd(eiv_comp_summary$SDM_estimate))
  
  eiv_comp_summary %>% 
    dplyr::select(-c(SDM_estimate, SDM_SE)) %>% 
    mutate(latent_lower = latent_estimate - 1.96*latent_SE,
           latent_upper = latent_estimate + 1.96*latent_SE) %>% 
    mutate(SDM_lower = SDM_estimate_scaled - 1.96*SDM_SE_scaled,
           SDM_upper = SDM_estimate_scaled + 1.96*SDM_SE_scaled) -> eiv_comp_summary
  
  return(eiv_comp_summary)
}

#### Plot EiV ####
compare_EiV_plot <- function(eiv_comp_summary, plot_title){
  plot <- ggplot(eiv_comp_summary, aes(x = SDM_estimate_scaled, y = latent_estimate, 
                                    xmin = latent_lower, xmax = latent_upper,
                                    ymin = SDM_lower, ymax = SDM_upper)) +
    geom_point() + 
    geom_errorbar(width = 0.1) +
    geom_errorbarh(aes(xmin = latent_lower, xmax = latent_upper), height = 0.1) +
    xlab("SDM Estimate (Data)") +
    ylab("SAR Estimate (Latent State)") +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed")
    
  
  return(plot)
}



#### Run for each model ####
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_csyif_index_t <- reformat_latent_state_df(SAR_summary = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_summary,
SDM_summary = seabird_SDM_summary,
SAR_years = common_years_seabird_prey,
SDM_years = jsoes_years,
SAR_parameter_name = "csyif_index_t_latent",
SDM_parameter_name = "csyif_index_of_abundance")

ggplot(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_csyif_index_t, aes(x = SDM_estimate_scaled, y = latent_estimate)) +
  geom_point()


#### EiV - Seabird x yearlings ####

## Step 1: Summarize + reformat the outputs from the SAR model and the SDM model for comparison
# SAR estimates + SE
load(here::here("R", "06_stage2_SAR", "06.2_seabird_SAR", "SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output_zscored.rda"))
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Obj <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Obj
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_output$SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt

as.data.frame(summary(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt$SD)) %>% 
  rownames_to_column("parameter") -> SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_summary

# SDM estimates + SE
load(here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "seabird_SDM_output.rda"))
seabird_SDM_Obj <- seabird_SDM_output$seabird_SDM_Obj
seabird_SDM_Opt <- seabird_SDM_output$seabird_SDM_Opt
seabird_SDM_report <- seabird_SDM_output$seabird_SDM_report

as.data.frame(summary(seabird_SDM_Opt$SD)) %>% 
  rownames_to_column("parameter") -> seabird_SDM_summary

subset(seabird_SDM_summary, grepl("csyif_index_of_abundance", parameter))



# csyif_index_t_latent
# pianka_o_csyif_sosh_t_latent
# sosh_index_t_latent
# pianka_o_csyif_comu_t_latent
# comu_index_t_latent
# pianka_o_csyif_prey_field_t_latent
# prey_field_index_t_latent



#### EiV - Seabird x subyearlings ####
# cssif_index_t_latent
# pianka_o_cssif_sosh_t_latent
# sosh_index_t_latent
# pianka_o_cssif_comu_t_latent
# comu_index_t_latent
# pianka_o_cssif_prey_field_t_latent
# prey_field_index_t_latent


#### EiV - hake x yearlings ####
# csyif_index_t_latent
# pianka_o_csyif_hake_t_latent
# hake_index_t_latent
# pianka_o_csyif_prey_field_t_latent
# prey_field_index_t_latent


#### EiV - hake x subyearlings ####
# cssif_index_t_latent
# pianka_o_cssif_hake_t_latent
# hake_index_t_latent
# pianka_o_cssif_prey_field_t_latent
# prey_field_index_t_latent






# CSYIF index
plot(x = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent, y = (SRF_v1_seabird_prey_predictors$csyif_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))

csyif_index_eiv_comp_data <- generate_eiv_comp_data(SDM_estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent/max(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent),
                                                    SDM_SE = sqrt(diag(csyif_aggregate_prey_field_seabird_model_Opt$SD$cov[csyif_index_of_abundance_cov_indices, csyif_index_of_abundance_cov_indices]/max(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_of_abundance)^2)),
                                                    SAR_estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent,
                                                    SAR_SE = subset(integrated_seabird_prey_field_SAR_SD_summary, grepl("csyif_index_t_latent", parameter))$std_error,
                                                    common_years = common_years_seabirds)

csyif_index_eiv_plot <- compare_EiV_plot(eiv_comp_data = csyif_index_eiv_comp_data, 
                                         covariate_name = "CSYIF June JSOES index of abundance")

ggsave(here::here("two_stage_models", "csyif_alternative_models_SWFSC", "seabird_model", "figures", "csyif_index_eiv_plot.png"),
       csyif_index_eiv_plot,
       height = 6, width = 8) 

