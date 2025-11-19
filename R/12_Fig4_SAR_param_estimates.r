## 12_Fig4_SAR_param_estimates

# Description: This script generates Figure 4, which shows the parameter estimates
# for the various SAR models


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


#### Compare fixed effects across models ####

effect_name_df <- data.frame(effect = c("hake_index", "ov_hake", "ov_prey_field", "prey_field_index",
                                        "sosh_index", "ov_sosh", "comu_index", "ov_comu"),
                             name = c("Hake\nAbundance", "Hake\nOverlap", "Prey\nOverlap", "Prey\nAbundance",
                                      "SOSH\nAbundance", "SOSH\nOverlap",
                                      "COMU\nAbundance", "COMU\nOverlap"))

effect_name_df$name <- factor(effect_name_df$name, 
                              levels = c("Prey\nAbundance", "Prey\nOverlap", 
                                         "Hake\nAbundance", "Hake\nOverlap", 
                                         "SOSH\nAbundance", "SOSH\nOverlap",
                                         "COMU\nAbundance", "COMU\nOverlap"))

# function to extract fixed effects
extract_fixed_effects_uncertainty <- function(model_opt_sd){
  as.data.frame(summary(model_opt_sd)) %>% 
    rownames_to_column("parameter") %>% 
    janitor::clean_names() -> SD_summary
  
  # extract the fixed effects
  fixed_effects_SD_summary <- subset(SD_summary, grepl("beta", parameter))
  
  fixed_effects_SD_summary %>% 
    mutate(upper = estimate + 1.96 * std_error,
           lower = estimate - 1.96 * std_error) %>% 
    mutate(significance = ifelse(lower > 0 | upper < 0, 
                                 "significant", "not significant")) %>% 
    mutate(effect = gsub("beta_", "", parameter)) %>% 
    left_join(effect_name_df, by = "effect") -> fixed_effects_SD_summary
  
  return(fixed_effects_SD_summary)
}

# function to visualize fixed effects estimates

plot_FE_estimates <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
                             drop_outmigration = TRUE,
                             drop_transport = TRUE,
                             drop_chinook_abundance = TRUE){
  if(drop_intercept == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, parameter != "beta_0")
  }
  if(drop_outmigration == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_outmigration", "beta_outmigration2")))
  }
  if(drop_transport == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, parameter != "beta_transport")
  }
  if(drop_chinook_abundance == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_csyif_index", "beta_cssif_index")))
  }
  
  significance_colors = c("significant" = "black",
                          "not significant" = "gray80")
  
  plot <- ggplot(fixed_effects_SD_summary, aes(x = name, y = estimate,
                                               ymax = upper, ymin = lower)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(size = 5) +
    geom_errorbar(width = 0.2) +
    # scale_color_manual(values = significance_colors) +
    guides(color = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.1, 0.1)) +
    xlab("Marine Survival Covariate") +
    ylab("Parameter Estimate") +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          legend.key.height = unit(1.25, "cm"),
          legend.key.width = unit(1.25, "cm"),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.title.x = element_text(size = 20, margin = margin(t = 10)),
          axis.title.y = element_text(size = 20, margin = margin(r = 10)),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
} 


### Extract the estimates from each of the model objects

## SRF models

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt$SD)
SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = SRF_06_3_hake_prey_csyif_only_SAR_SEcov_Opt$SD)
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_Opt$SD)
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_Opt$SD)

## UCSF models
UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE <- extract_fixed_effects_uncertainty(model_opt_sd = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_Opt$SD)
UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE <- extract_fixed_effects_uncertainty(model_opt_sd = UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_Opt$SD)
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE <- extract_fixed_effects_uncertainty(model_opt_sd = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_Opt$SD)
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE <- extract_fixed_effects_uncertainty(model_opt_sd = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_Opt$SD)


### Visualize the parameter estimates

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE_plot <- plot_FE_estimates(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE)
SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE_plot <- plot_FE_estimates(SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE)
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE_plot <- plot_FE_estimates(SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE)
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE_plot <- plot_FE_estimates(SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE)

UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE_plot <- plot_FE_estimates(UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE)
UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE_plot <- plot_FE_estimates(UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE)
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE_plot <- plot_FE_estimates(UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE)
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE_plot <- plot_FE_estimates(UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE)


seabird_FE_comparison_plot <- ggarrange(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE_plot, 
                                        SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE_plot,
                                        UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE_plot,
                                        UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE_plot,
                                        labels = c("SRF CSYIF", "SRF CSSIF", "UCSF CSYIF", "UCSF CSSIF"))

ggsave(here::here("figures", "paper_figures", "parameter_estimate_plots", "seabird_FE_comparison_plot.png"), seabird_FE_comparison_plot,  
       height = 12, width = 16)

hake_FE_comparison_plot <- ggarrange(SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE_plot, 
                                     SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE_plot,
                                     UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE_plot,
                                     UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE_plot,
                                        labels = c("SRF CSYIF", "SRF CSSIF", "UCSF CSYIF", "UCSF CSSIF"))

ggsave(here::here("figures", "paper_figures", "parameter_estimate_plots", "hake_FE_comparison_plot.png"), hake_FE_comparison_plot,  
       height = 12, width = 16)


#### Plot model fit to data ####

# function to sample from estimate + SE
generate_samples <- function(estimate, std_error, nsamples = 1000){
  param_samples <- rnorm(nsamples, mean = estimate, sd = std_error)
  return(param_samples)
}

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples <- data.frame(
  beta_0 = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_0","estimate"],
                   std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_0","std_error"]),
  beta_csyif_index = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"],
                                      std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index","std_error"]),
  
  beta_prey_field_index = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"],
                            std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index","std_error"]),
  beta_ov_prey_field = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"],
                            std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field","std_error"]),
  beta_sosh_index = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"],
                            std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_sosh_index","std_error"]),
  beta_ov_sosh = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"],
                            std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh","std_error"]),
  beta_comu_index = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"],
                            std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_comu_index","std_error"]),
  beta_ov_comu = generate_samples(estimate = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"],
                            std_error = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_comu","std_error"])
  )

# get predictors in a df
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                       csyif_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent,
                                                                       pianka_o_csyif_sosh_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_sosh_t_latent,
                                                                       sosh_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$sosh_index_t_latent,
                                                                       pianka_o_csyif_comu_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_comu_t_latent,
                                                                       comu_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$comu_index_t_latent,
                                                                       pianka_o_csyif_prey_field_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_prey_field_t_latent,
                                                                       prey_field_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$prey_field_index_t_latent)


SRF_06_2_seabird_prey_csyif_only_SAR_prob_samples <- as.data.frame(matrix(nrow = 1000, ncol = nrow(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df)))
colnames(SRF_06_2_seabird_prey_csyif_only_SAR_prob_samples) <- common_years_seabird_prey

# Loop through samples to generate distribution of predicted probability
for (i in 1:nrow(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df)){ # loop through each of the years
  for (j in 1:nrow(SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples)){ # loop through each of the samples
    
    linear_predictor <- 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_0"] +
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_csyif_index"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$csyif_index_t_latent[i] +
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_prey_field_index"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent[i] +
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_ov_prey_field"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent[i] +
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_sosh_index"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$sosh_index_t_latent[i] +
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_ov_sosh"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_sosh_t_latent[i] +
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_comu_index"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$comu_index_t_latent[i] +
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_param_samples[j, "beta_ov_comu"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_comu_t_latent[i]
   
    SRF_06_2_seabird_prey_csyif_only_SAR_prob_samples[j,i] <- inv.logit(linear_predictor) 
  }
}

# get median + 95% confidence interval from these samples
apply(SRF_06_2_seabird_prey_csyif_only_SAR_prob_samples, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
  as.data.frame() %>% 
  rownames_to_column("quantile") %>% 
  pivot_longer(cols = -c("quantile"), names_to = "run_year", values_to = "prob") %>% 
  mutate(run_year = as.numeric(run_year)) -> SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles

SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles %>% 
  pivot_wider(names_from = "quantile", values_from = "prob") -> SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles

## SRF CSYIF seabird as an example model to visualize

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_Opt$SD
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report

# plot MLE
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                       csyif_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent,
                            pianka_o_csyif_sosh_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_sosh_t_latent,
                            sosh_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$sosh_index_t_latent,
                            pianka_o_csyif_comu_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_comu_t_latent,
                            comu_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$comu_index_t_latent,
                            pianka_o_csyif_prey_field_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_prey_field_t_latent,
                            prey_field_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$prey_field_index_t_latent)

# beta_0 + 
#   // individual level covariates
# beta_transport * transport_i(i) +
#   beta_outmigration * outmigration_i(i) +
#   beta_outmigration2 * pow(outmigration_i(i),2) +
#   // SDM derived outputs
# beta_csyif_index*csyif_index_t_latent(t_i(i)) +
#   beta_prey_field_index*prey_field_index_t_latent(t_i(i)) +
#   beta_ov_prey_field*pianka_o_csyif_prey_field_t_latent(t_i(i)) +
#   beta_sosh_index*sosh_index_t_latent(t_i(i)) +
#   beta_ov_sosh*pianka_o_csyif_sosh_t_latent(t_i(i)) +
#   beta_comu_index*comu_index_t_latent(t_i(i)) +
#   beta_ov_comu*pianka_o_csyif_comu_t_latent(t_i(i)),


SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
  linear_predictor = 
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_0","estimate"] +
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
      SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$csyif_index_t_latent +
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$sosh_index_t_latent +
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_sosh_t_latent +
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$comu_index_t_latent +
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
    SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_comu_t_latent)

SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions$linear_predictor)

SRF_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_subset_SAR


SRF_subset_SAR %>% 
  left_join(SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions, by = "run_year") -> SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions_comp


SRF_06_2_seabird_prey_csyif_only_SAR_MLE_predict_plot <- ggplot(SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
  scale_shape_manual(name = NULL, values = c("Empirical" = 19, "Predicted" = 8)) +
  xlab("Run Year") +
  ylab("Marine Survival") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))

ggsave(here::here("figures", "paper_figures", "SRF_06_2_seabird_prey_csyif_only_SAR_MLE_predict_plot.png"), SRF_06_2_seabird_prey_csyif_only_SAR_MLE_predict_plot,  
       height = 6, width = 8)


# Show uncertainty from estimation using sampling
SRF_06_2_seabird_prey_csyif_only_SAR_predict_plot <- ggplot(SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`, ymax = `97.5%`, ymin = `2.5%`), width = 0.1) +
  # geom_ribbon(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`, ymax = `97.5%`, ymin = `2.5%`), fill = "grey70") +
  # geom_linerange(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`, ymax = `97.5%`, ymin = `2.5%`), fill = "grey70", linetype = 2) +
  geom_point() +
  geom_point(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`), shape = 8)

SRF_06_2_seabird_prey_csyif_only_SAR_predict_plot <- ggplot(SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`, ymax = `97.5%`, ymin = `2.5%`), width = 0.1) +
  # geom_ribbon(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`, ymax = `97.5%`, ymin = `2.5%`), fill = "grey70") +
  # geom_linerange(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`, ymax = `97.5%`, ymin = `2.5%`), fill = "grey70", linetype = 2) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = SRF_06_2_seabird_prey_csyif_only_SAR_prob_quantiles, aes(x = run_year, y = `50%`, shape = "Predicted")) +
  scale_shape_manual(name = NULL, values = c("Empirical" = 19, "Predicted" = 8)) +
  scale_y_continuous(breaks = seq(0, 0.35, 0.05)) +
  xlab("Run Year") +
  ylab("Marine Survival") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))

ggsave(here::here("figures", "paper_figures", "SRF_06_2_seabird_prey_csyif_only_SAR_predict_plot.png"), SRF_06_2_seabird_prey_csyif_only_SAR_predict_plot,  
       height = 6, width = 8)

