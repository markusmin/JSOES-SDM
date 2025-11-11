#### 07_model_output_comparison_figure

# This script will compare the parameter estimates from multiple models.

library(tidyverse)
library(here)
library(TMB)
library(ggforce)
library(sdmTMB)
library(kableExtra)
library(fmesher)
library(ggrepel)
# knitr::opts_chunk()

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

#### Load the Snake River Fall models ####

## 06.1 prey model

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.1_prey_field_SAR", "SRF_06_1_prey_field_SAR_csyif_only_output_zscored.rda"))

SRF_06_1_prey_field_SAR_csyif_only_Obj <- SRF_06_1_prey_field_SAR_csyif_only_output$SRF_06_1_prey_field_SAR_csyif_only_Obj
SRF_06_1_prey_field_SAR_csyif_only_Opt <- SRF_06_1_prey_field_SAR_csyif_only_output$SRF_06_1_prey_field_SAR_csyif_only_Opt
SRF_06_1_prey_field_SAR_csyif_only_report <- SRF_06_1_prey_field_SAR_csyif_only_output$SRF_06_1_prey_field_SAR_csyif_only_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.1_prey_field_SAR", "SRF_06_1_prey_field_SAR_cssif_only_output_zscored.rda"))

SRF_06_1_prey_field_SAR_cssif_only_Obj <- SRF_06_1_prey_field_SAR_cssif_only_output$SRF_06_1_prey_field_SAR_cssif_only_Obj
SRF_06_1_prey_field_SAR_cssif_only_Opt <- SRF_06_1_prey_field_SAR_cssif_only_output$SRF_06_1_prey_field_SAR_cssif_only_Opt
SRF_06_1_prey_field_SAR_cssif_only_report <- SRF_06_1_prey_field_SAR_cssif_only_output$SRF_06_1_prey_field_SAR_cssif_only_report

## 06.2 seabird model

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.2_seabird_SAR", "SRF_06_2_seabird_prey_csyif_only_SAR_output_zscored.rda"))

SRF_06_2_seabird_prey_csyif_only_SAR_Obj <- SRF_06_2_seabird_prey_csyif_only_SAR_output$SRF_06_2_seabird_prey_csyif_only_SAR_Obj
SRF_06_2_seabird_prey_csyif_only_SAR_Opt <- SRF_06_2_seabird_prey_csyif_only_SAR_output$SRF_06_2_seabird_prey_csyif_only_SAR_Opt
SRF_06_2_seabird_prey_csyif_only_SAR_report <- SRF_06_2_seabird_prey_csyif_only_SAR_output$SRF_06_2_seabird_prey_csyif_only_SAR_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.2_seabird_SAR", "SRF_06_2_seabird_prey_cssif_only_SAR_output_zscored.rda"))

SRF_06_2_seabird_prey_cssif_only_SAR_Obj <- SRF_06_2_seabird_prey_cssif_only_SAR_output$SRF_06_2_seabird_prey_cssif_only_SAR_Obj
SRF_06_2_seabird_prey_cssif_only_SAR_Opt = SRF_06_2_seabird_prey_cssif_only_SAR_output$SRF_06_2_seabird_prey_cssif_only_SAR_Opt
SRF_06_2_seabird_prey_cssif_only_SAR_report = SRF_06_2_seabird_prey_cssif_only_SAR_output$SRF_06_2_seabird_prey_cssif_only_SAR_report



## 06.3 hake model

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.3_hake_SAR", "SRF_06_3_hake_prey_csyif_only_SAR_output_zscored.rda"))

SRF_06_3_hake_prey_csyif_only_SAR_Obj <- SRF_06_3_hake_prey_csyif_only_SAR_output$SRF_06_3_hake_prey_csyif_only_SAR_Obj
SRF_06_3_hake_prey_csyif_only_SAR_Opt <- SRF_06_3_hake_prey_csyif_only_SAR_output$SRF_06_3_hake_prey_csyif_only_SAR_Opt
SRF_06_3_hake_prey_csyif_only_SAR_report <- SRF_06_3_hake_prey_csyif_only_SAR_output$SRF_06_3_hake_prey_csyif_only_SAR_report

# CSSIF

## NOTE: This model did not converge
load(here::here("R", "06_stage2_SAR", "06.3_hake_SAR", "SRF_06_3_hake_prey_cssif_only_SAR_output_zscored.rda"))

SRF_06_3_hake_prey_cssif_only_SAR_Obj <- SRF_06_3_hake_prey_cssif_only_SAR_output$SRF_06_3_hake_prey_cssif_only_SAR_Obj
SRF_06_3_hake_prey_cssif_only_SAR_Opt <- SRF_06_3_hake_prey_cssif_only_SAR_output$SRF_06_3_hake_prey_cssif_only_SAR_Opt
SRF_06_3_hake_prey_cssif_only_SAR_report <- SRF_06_3_hake_prey_cssif_only_SAR_output$SRF_06_3_hake_prey_cssif_only_SAR_report


#### Examine model estimates ####

SRF_06_1_prey_field_SAR_csyif_only_Opt$SD
SRF_06_1_prey_field_SAR_cssif_only_Opt$SD
SRF_06_2_seabird_prey_csyif_only_SAR_Opt$SD
SRF_06_2_seabird_prey_cssif_only_SAR_Opt$SD
SRF_06_3_hake_prey_csyif_only_SAR_Opt$SD
SRF_06_3_hake_prey_cssif_only_SAR_Opt$SD





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

plot_FE_estmates <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
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


SRF_06_1_prey_field_SAR_csyif_only_FE <- extract_fixed_effects_uncertainty(model_opt_sd = SRF_06_1_prey_field_SAR_csyif_only_Opt$SD)

SRF_06_2_seabird_prey_csyif_only_SAR_FE <- extract_fixed_effects_uncertainty(model_opt_sd = SRF_06_2_seabird_prey_csyif_only_SAR_Opt$SD)

SRF_06_3_hake_prey_csyif_only_SAR_FE <- extract_fixed_effects_uncertainty(model_opt_sd = SRF_06_3_hake_prey_csyif_only_SAR_Opt$SD)


### Visualize the parameter estimates

prey_field_FE_plot <- plot_FE_estmates(SRF_06_1_prey_field_SAR_csyif_only_FE)

seabird_FE_plot <- plot_FE_estmates(SRF_06_2_seabird_prey_csyif_only_SAR_FE)

hake_FE_plot <- plot_FE_estmates(SRF_06_3_hake_prey_csyif_only_SAR_FE)

ggsave(here::here("figures", "presentation_figures", "seabird_FE_plot.png"), seabird_FE_plot,  
       height = 6, width = 8)

ggsave(here::here("figures", "presentation_figures", "hake_FE_plot.png"), hake_FE_plot,  
       height = 6, width = 8)


### combine the model outputs together

SRF_06_1_prey_field_SAR_csyif_only_subset <- subset(fixed_effects_SD_summary_SRF_06_1_prey_field_SAR_csyif_only, parameter %in% c("beta_csyif_index", "beta_prey_field_index","beta_ov_prey_field"))

seabird_model_subset <- subset(fixed_effects_SD_summary_seabird_model, parameter %in% c(
                                                                "beta_sosh_index", "beta_comu_index",
                                                                 "beta_ov_sosh", "beta_ov_comu"))

hake_model_subset <- subset(fixed_effects_SD_summary_hake_model, parameter %in% c("beta_hake_index", "beta_ov_hake"))

SRF_06_1_prey_field_SAR_csyif_only_subset %>% 
  bind_rows(seabird_model_subset) %>% 
  bind_rows(hake_model_subset) -> model_fixed_effects_param_comparison

model_fixed_effects_param_comparison %>% 
  mutate(effect = gsub("beta_", "", parameter)) -> model_fixed_effects_param_comparison

model_fixed_effects_param_comparison %>% 
  mutate(significance = ifelse(estimate - 1.96 * std_error > 0 |
                                estimate + 1.96 < 0, 
                              "significant", "not significant")) -> model_fixed_effects_param_comparison

significance_colors = c("significant" = "black",
                       "not significant" = "gray80")

model_fixed_effects_param_comparison$effect <- factor(model_fixed_effects_param_comparison$effect, 
                                                         levels = c("csyif_index",
                                                                    "prey_field_index",
                                                                    "ov_prey_field",
                                                                    "comu_index",
                                                                    "ov_comu",
                                                                    "sosh_index",
                                                                    "ov_sosh",
                                                                    "hake_index",
                                                                    "ov_hake"
                                                                    ))

# drop the csyif index for plotting, since we're focusing on ecological interactions
model_fixed_effects_param_comparison <- subset(model_fixed_effects_param_comparison, effect != "csyif_index")

model_fixed_effects_param_comparison_plot <- ggplot(model_fixed_effects_param_comparison, aes(x = effect, y = estimate, color = significance,
                                                 ymax = estimate + 1.96*std_error, ymin = estimate - 1.96 * std_error)) +
  geom_point() +
  geom_errorbar() +
  scale_color_manual(values = significance_colors) +
  guides(color = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.1, 0.1))

ggsave(here::here("figures", "SAR_model", "model_fixed_effects_param_comparison_plot.png"),
       model_fixed_effects_param_comparison_plot,
       height = 6, width = 10) 


#### Visualize the prey field results ####

# use the standard errors to get parameters with uncertainty
as.data.frame(summary(SRF_06_1_prey_field_SAR_v1_Opt$SD)) %>% 
  rownames_to_column("parameter") %>% 
  janitor::clean_names() -> SRF_06_1_prey_field_SAR_v1_SD_summary

# extract the fixed effects
fixed_effects_SD_summary_prey_field_model <- subset(SRF_06_1_prey_field_SAR_v1_SD_summary, grepl("beta", parameter))

# extract the values of the latent variables from the EiV approach
latent_variables <- c("csyif_index_t_latent", "cssif_index_t_latent",
                      "prey_field_index_t_latent", "pianka_o_csyif_prey_field_t_latent",
                      "pianka_o_csyif_prey_field_t_latent")

csyif_index_t_latent <- subset(SRF_06_1_prey_field_SAR_v1_SD_summary, grepl("csyif_index_t_latent", parameter))
cssif_index_t_latent <- subset(SRF_06_1_prey_field_SAR_v1_SD_summary, grepl("cssif_index_t_latent", parameter))
prey_field_index_t_latent <- subset(SRF_06_1_prey_field_SAR_v1_SD_summary, grepl("prey_field_index_t_latent", parameter))
pianka_o_csyif_prey_field_t_latent <- subset(SRF_06_1_prey_field_SAR_v1_SD_summary, grepl("pianka_o_csyif_prey_field_t_latent", parameter))
pianka_o_cssif_prey_field_t_latent <- subset(SRF_06_1_prey_field_SAR_v1_SD_summary, grepl("pianka_o_cssif_prey_field_t_latent", parameter))

# are the overlap values correlated for csyif and cssif?
# Oh damn, super
plot(x = pianka_o_csyif_prey_field_t_latent$estimate,
     y = pianka_o_cssif_prey_field_t_latent$estimate)

cor(pianka_o_csyif_prey_field_t_latent$estimate, 
    pianka_o_cssif_prey_field_t_latent$estimate)


# add the mean covariate values for all fixed effects
fixed_effects_SD_summary_prey_field_model$mean_covariate_value <- c(1, 
                                                                 0, 
                                                                 0, 
                                                                 0,
                                                                 mean(csyif_index_t_latent$estimate),
                                                                 mean(cssif_index_t_latent$estimate),
                                                                 mean(prey_field_index_t_latent$estimate),
                                                                 mean(pianka_o_csyif_prey_field_t_latent$estimate),
                                                                 mean(pianka_o_cssif_prey_field_t_latent$estimate))

fixed_effects_SD_summary_prey_field_model %>% 
  mutate(mean_effect = estimate * mean_covariate_value) -> fixed_effects_SD_summary_prey_field_model


# CSYIF Index
beta_csyif_index_post <- generate_posterior_predictive(SD_summary = SRF_06_1_prey_field_SAR_v1_SD_summary,
                                                       fixed_effects_SD_summary = fixed_effects_SD_summary_prey_field_model,
                                                       parameter_name = "beta_csyif_index", covariate_values = seq(0, 1, by = 0.01))
csyif_index_post_pred_plot <- plot_posterior_predictive(beta_csyif_index_post, "CSYIF June JSOES index of abundance")
ggsave(here::here("figures", "SAR_model", "prey_field_model", "csyif_index_post_pred_plot.png"),
       csyif_index_post_pred_plot,
       height = 6, width = 8) 

# CSSIF Index
beta_cssif_index_post <- generate_posterior_predictive(SD_summary = SRF_06_1_prey_field_SAR_v1_SD_summary,
                                                       fixed_effects_SD_summary = fixed_effects_SD_summary_prey_field_model,
                                                       parameter_name = "beta_cssif_index", covariate_values = seq(0, 1, by = 0.01))
cssif_index_post_pred_plot <- plot_posterior_predictive(beta_cssif_index_post, "CSSIF June JSOES index of abundance")
ggsave(here::here("figures", "SAR_model", "prey_field_model", "cssif_index_post_pred_plot.png"),
       cssif_index_post_pred_plot,
       height = 6, width = 8) 

# prey_field_index
beta_prey_field_index_post <- generate_posterior_predictive(SD_summary = SRF_06_1_prey_field_SAR_v1_SD_summary,
                                                            fixed_effects_SD_summary = fixed_effects_SD_summary_prey_field_model,
                                                            parameter_name = "beta_prey_field_index", covariate_values = seq(0, 1, by = 0.01))
prey_field_index_post_pred_plot <- plot_posterior_predictive(beta_prey_field_index_post, "Aggregate prey field index of abundance")
ggsave(here::here("figures", "SAR_model", "prey_field_model", "prey_field_index_post_pred_plot.png"),
       prey_field_index_post_pred_plot,
       height = 6, width = 8) 


# ov_csyif_prey_field
beta_ov_csyif_prey_field_post <- generate_posterior_predictive(SD_summary = SRF_06_1_prey_field_SAR_v1_SD_summary,
                                                         fixed_effects_SD_summary = fixed_effects_SD_summary_prey_field_model,
                                                         parameter_name = "beta_ov_csyif_prey_field", covariate_values = seq(0, 1, by = 0.01))
ov_csyif_prey_field_post_pred_plot <- plot_posterior_predictive(beta_ov_csyif_prey_field_post, "CSYIF x prey field overlap")
ggsave(here::here("figures", "SAR_model", "prey_field_model", "ov_csyif_prey_field_post_pred_plot.png"),
       ov_csyif_prey_field_post_pred_plot,
       height = 6, width = 8) 

# ov_cssif_prey_field
beta_ov_cssif_prey_field_post <- generate_posterior_predictive(SD_summary = SRF_06_1_prey_field_SAR_v1_SD_summary,
                                                               fixed_effects_SD_summary = fixed_effects_SD_summary_prey_field_model,
                                                               parameter_name = "beta_ov_cssif_prey_field", covariate_values = seq(0, 1, by = 0.01))
ov_cssif_prey_field_post_pred_plot <- plot_posterior_predictive(beta_ov_cssif_prey_field_post, "CSSIF x prey field overlap")
ggsave(here::here("figures", "SAR_model", "prey_field_model", "ov_cssif_prey_field_post_pred_plot.png"),
       ov_cssif_prey_field_post_pred_plot,
       height = 6, width = 8) 



outmigration_post <- generate_posterior_predictive_outmigration(SD_summary = SRF_06_1_prey_field_SAR_v1_SD_summary,
                                                                fixed_effects_SD_summary = fixed_effects_SD_summary_prey_field_model,
                                                                covariate_values = seq(min(SRF_06_1_prey_field_SAR_v1_Data$outmigration_i), max(SRF_06_1_prey_field_SAR_v1_Data$outmigration_i), length.out = 100))
outmigration_post_pred_plot <- plot_posterior_predictive(outmigration_post, "Date of outmigration")

ggsave(here::here("figures", "SAR_model", "prey_field_model", "outmigration_post_pred_plot.png"),
       outmigration_post_pred_plot,
       height = 6, width = 8) 







