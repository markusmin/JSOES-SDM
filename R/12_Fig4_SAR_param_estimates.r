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


#### Version for paper figure

plot_FE_estimates_forfig <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
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
    facet_wrap(~model, ncol = 1) +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          legend.key.height = unit(1.25, "cm"),
          legend.key.width = unit(1.25, "cm"),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 12),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          axis.title.x = element_text(size = 20, margin = margin(t = 10)),
          axis.title.y = element_text(size = 20, margin = margin(r = 10)),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
} 

plot_FE_estimates_forfig_combined <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
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
    facet_wrap(~model, ncol = 2, scales = "free_x") +
    # facet_wrap(~model, ncol = 2) +
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

# SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$model <- "Snake River Fall - Yearlings x Seabirds"
# SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$model <- "Snake River Fall - Yearlings x Hake"
# SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$model  <- "Snake River Fall - Subyearlings x Seabirds"
# SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$model  <- "Snake River Fall - Subyearlings x Hake"
# 
# UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Yearlings x Seabirds"
# UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Yearlings x Hake"
# UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Subyearlings x Seabirds"
# UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Subyearlings x Hake"

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$model <- "Snake River Fall - Yearlings"
SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$model <- "Snake River Fall - Yearlings"
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$model  <- "Snake River Fall - Subyearlings"
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$model  <- "Snake River Fall - Subyearlings"

UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Yearlings"
UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Yearlings"
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Subyearlings"
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$model <- "Upper Columbia Summer/Fall - Subyearlings"


SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE %>% 
  bind_rows(SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE) %>% 
  bind_rows(SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE) %>% 
  bind_rows(SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE) %>% 
  bind_rows(UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE) %>% 
  bind_rows(UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE) %>% 
  bind_rows(UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE) %>% 
  bind_rows(UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE) -> SAR_models_FE

SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE %>% 
  bind_rows(SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE) %>% 
  bind_rows(UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE) %>% 
  bind_rows(UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE) -> SAR_hake_models_FE

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE %>% 
  bind_rows(SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE) %>% 
  bind_rows(UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE) %>% 
  bind_rows(UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE) -> SAR_seabird_models_FE

SAR_seabird_models_FE_plot <- plot_FE_estimates_forfig(SAR_seabird_models_FE)
SAR_hake_models_FE_plot <- plot_FE_estimates_forfig(SAR_hake_models_FE)

ggsave(here::here("figures", "paper_figures", "parameter_estimate_plots", "SAR_seabird_models_FE_plot.png"), SAR_seabird_models_FE_plot,  
       height = 10, width = 6)

ggsave(here::here("figures", "paper_figures", "parameter_estimate_plots", "SAR_hake_models_FE_plot.png"), SAR_hake_models_FE_plot,  
       height = 10, width = 6)

fig4_SAR_parameter_estimates <- ggarrange(SAR_seabird_models_FE_plot,
                                          SAR_hake_models_FE_plot,
                                   labels = c("(A)", "(B)"),
                                   font.label = list(size = 20, face = "plain"),
                                   label.x = 0, label.y = 1,
                                   widths = c(1.2, 1),
                                   ncol = 2)

ggsave(here::here("figures", "paper_figures", "fig4_SAR_parameter_estimates.png"), fig4_SAR_parameter_estimates,  
       height = 12, width = 12)

SAR_FE_plot <- plot_FE_estimates_forfig_combined(SAR_models_FE)

ggsave(here::here("figures", "paper_figures", "parameter_estimate_plots", "fig4_SAR_parameter_estimates_v2.png"), SAR_FE_plot,  
       height = 12, width = 12)


### New version with labels for each panel


#### Plot model fit to data ####

# prep SAR data

#### Snake River Fall Chinook

# extract Snake River Fall Chinook from SAR data
SRF <- subset(SAR_data, run_name == "Fall" & group == "Snake River")

#### subset by outmigration date
hist(SRF$outmigration_date)
summary(SRF$outmigration_date)

# restrict the analysis to fish that could have theoretically been caught in the June JSOES survey.
# we will call this April 15 - June 15 (for now!)
# day 105 is April 15 in a non-leap year and day 167 is June 15 in a leap year
# SRF_subset <- filter(SRF, outmigration_date >= 105 & outmigration_date <= 167)

# let's just have a final cutoff - based on the fact that fish can hang out
# off the coast for a while, it's less justified to have a front-end cutoff date
SRF_subset <- filter(SRF, outmigration_date <= 167)
nrow(SRF_subset)/nrow(SRF)
# this leaves us with 54% of our initial dataset


#### Join with the transport data

transport_data_export <- read.csv(here::here("model_inputs", "chinook_transport.csv"))

# join the juvenile transport data with the SAR data
SRF_subset %>% 
  left_join(transport_data_export, join_by(tag_code == tag_id)) -> SRF_subset

SRF_subset$transport[is.na(SRF_subset$transport)] <- 0

# inspect number of transported fish
table(SRF_subset$transport)


#### Drop the few natural origin fish

# fish with rear type U or H are hatchery fish
SRF_subset %>% 
  mutate(rear_numeric = ifelse(rear_type_code %in% c("H", "U"), 1, 0)) -> SRF_subset

# ok, so they're all hatchery. That's perhaps not surprising, given what we see in the JSOES survey is basically all hatchery
# Let's only keep hatchery fish
SRF_subset <- subset(SRF_subset, rear_type_code != "W")


#### Upper Columbia Summer/Fall Chinook

# extract Snake River Fall Chinook from SAR data
UCSF <- subset(SAR_data, run_name %in% c("Summer", "Fall") & group == "Upper Columbia")

#### subset by outmigration date
hist(UCSF$outmigration_date)
summary(UCSF$outmigration_date)

# restrict the analysis to fish that could have theoretically been caught in the June JSOES survey.
# we will call this April 15 - June 15 (for now!)
# day 105 is April 15 in a non-leap year and day 167 is June 15 in a leap year
# UCSF_subset <- filter(UCSF, outmigration_date >= 105 & outmigration_date <= 167)

# let's just have a final cutoff - based on the fact that fish can hang out
# off the coast for a while, it's less justified to have a front-end cutoff date
UCSF_subset <- filter(UCSF, outmigration_date <= 167)
nrow(UCSF_subset)/nrow(UCSF)
# this leaves us with 72% of our initial dataset

#### Drop the few natural origin fish (28 of ~97,500)

table(UCSF_subset$rear_type_code)

# fish with rear type U or H are hatchery fish
UCSF_subset %>% 
  mutate(rear_numeric = ifelse(rear_type_code %in% c("H", "U"), 1, 0)) -> UCSF_subset

# ok, so they're all hatchery. That's perhaps not surprising, given what we see in the JSOES survey is basically all hatchery
# Let's only keep hatchery fish
UCSF_subset <- subset(UCSF_subset, rear_type_code != "W")



#### function for analytical variance calculation
var_2_rv <- function(est1, var1, est2, var2){
  total_var <- var1*var2 + var1*est2^2 + var2*est1^2
  return(total_var)
}


# look at our latent estimates vs. the input data
# ok good, it's recovering the input well
plot(x = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent, y = (SRF_v1_seabird_prey_predictors$csyif_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))


#### SRF - Yearlings x Seabirds ####
SRF_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_subset_SAR

SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$variance <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$std_error^2

## variance of beta_0 parameter
var_beta_0 <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_0", "variance"]

## csyif index of abundance
# expected value of beta_csyif index parameter
est_beta_csyif_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index", "estimate"]

# variance of beta_csyif index parameter
var_beta_csyif_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index", "variance"]

# expected value of csyif_index of abundance predictor
est_csyif_index_of_abundance <- (SRF_v1_seabird_prey_predictors$csyif_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance)

# variance of csyif index of abundance predictor
var_csyif_index_of_abundance <- sqrt(diag(csyif_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (SRF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t - mean(SRF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t))/sd(SRF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_csyif_prey_field_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t))^2))


## sosh index
# expected value of beta_sosh index parameter
est_beta_sosh_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_sosh_index", "estimate"]

# variance of beta_sosh index parameter
var_beta_sosh_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_sosh_index", "variance"]

# expected value of sosh_index of abundance predictor
est_sosh_index_of_abundance <- (SRF_v1_seabird_prey_predictors$sosh_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$sosh_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$sosh_index_of_abundance)

# variance of sosh index of abundance predictor
var_sosh_index_of_abundance <- sqrt(diag(sosh_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$sosh_index_of_abundance))^2))

## sosh overlap
# expected value of beta_sosh parameter
est_beta_ov_sosh <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh", "estimate"]

# variance of beta_sosh parameter
var_beta_ov_sosh <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh", "variance"]

# expected value of sosh overlap predictor
est_sosh_overlap <- (SRF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t - mean(SRF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t))/sd(SRF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t)

# variance of sosh index of abundance predictor
var_sosh_overlap <- sqrt(diag(pianka_o_csyif_sosh_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t))^2))

## comu index
# expected value of beta_comu index parameter
est_beta_comu_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_comu_index", "estimate"]

# variance of beta_comu index parameter
var_beta_comu_index <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_comu_index", "variance"]

# expected value of comu_index of abundance predictor
est_comu_index_of_abundance <- (SRF_v1_seabird_prey_predictors$comu_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$comu_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$comu_index_of_abundance)

# variance of comu index of abundance predictor
var_comu_index_of_abundance <- sqrt(diag(comu_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$comu_index_of_abundance))^2))

## comu overlap
# expected value of beta_comu parameter
est_beta_ov_comu <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_comu", "estimate"]

# variance of beta_comu parameter
var_beta_ov_comu <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_comu", "variance"]

# expected value of comu overlap predictor
est_comu_overlap <- (SRF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t - mean(SRF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t))/sd(SRF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t)

# variance of comu index of abundance predictor
var_comu_overlap <- sqrt(diag(pianka_o_csyif_comu_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t))^2))


#### Calculate and plot total variance by year
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
           est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                       csyif_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent,
                                                                       pianka_o_csyif_sosh_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_sosh_t_latent,
                                                                       sosh_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$sosh_index_t_latent,
                                                                       pianka_o_csyif_comu_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_comu_t_latent,
                                                                       comu_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$comu_index_t_latent,
                                                                       pianka_o_csyif_prey_field_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_prey_field_t_latent,
                                                                       prey_field_index_t_latent = SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_report$prey_field_index_t_latent)

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

SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions$var <- SRF_06_2_seabird_prey_csyif_only_SAR_SEcov_pred_var

SRF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions

SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions %>% 
  left_join(SRF_subset_SAR, by = "run_year") -> SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp


SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp_plot <- ggplot(SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp_plot.png"), SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### SRF - Subyearlings x Seabirds ####
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$variance <-SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$std_error^2

## variance of beta_0 parameter
var_beta_0 <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_0", "variance"]

## cssif index of abundance
# expected value of beta_cssif index parameter
est_beta_cssif_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_cssif_index", "estimate"]

# variance of beta_cssif index parameter
var_beta_cssif_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_cssif_index", "variance"]

# expected value of cssif_index of abundance predictor
est_cssif_index_of_abundance <- (SRF_v1_seabird_prey_predictors$cssif_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$cssif_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$cssif_index_of_abundance)

# variance of cssif index of abundance predictor
var_cssif_index_of_abundance <- sqrt(diag(cssif_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$cssif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (SRF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t - mean(SRF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t))/sd(SRF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_cssif_prey_field_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t))^2))


## sosh index
# expected value of beta_sosh index parameter
est_beta_sosh_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_sosh_index", "estimate"]

# variance of beta_sosh index parameter
var_beta_sosh_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_sosh_index", "variance"]

# expected value of sosh_index of abundance predictor
est_sosh_index_of_abundance <- (SRF_v1_seabird_prey_predictors$sosh_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$sosh_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$sosh_index_of_abundance)

# variance of sosh index of abundance predictor
var_sosh_index_of_abundance <- sqrt(diag(sosh_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$sosh_index_of_abundance))^2))

## sosh overlap
# expected value of beta_sosh parameter
est_beta_ov_sosh <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh", "estimate"]

# variance of beta_sosh parameter
var_beta_ov_sosh <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh", "variance"]

# expected value of sosh overlap predictor
est_sosh_overlap <- (SRF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t - mean(SRF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t))/sd(SRF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t)

# variance of sosh index of abundance predictor
var_sosh_overlap <- sqrt(diag(pianka_o_cssif_sosh_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t))^2))

## comu index
# expected value of beta_comu index parameter
est_beta_comu_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_comu_index", "estimate"]

# variance of beta_comu index parameter
var_beta_comu_index <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_comu_index", "variance"]

# expected value of comu_index of abundance predictor
est_comu_index_of_abundance <- (SRF_v1_seabird_prey_predictors$comu_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$comu_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$comu_index_of_abundance)

# variance of comu index of abundance predictor
var_comu_index_of_abundance <- sqrt(diag(comu_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$comu_index_of_abundance))^2))

## comu overlap
# expected value of beta_comu parameter
est_beta_ov_comu <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_comu", "estimate"]

# variance of beta_comu parameter
var_beta_ov_comu <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_comu", "variance"]

# expected value of comu overlap predictor
est_comu_overlap <- (SRF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t - mean(SRF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t))/sd(SRF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t)

# variance of comu index of abundance predictor
var_comu_overlap <- sqrt(diag(pianka_o_cssif_comu_cov_matrix_common_years_seabird_prey * (1/sd(SRF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t))^2))

#### Calculate and plot total variance by year
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
           est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                       cssif_index_t_latent = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report$cssif_index_t_latent,
                                                                       pianka_o_cssif_sosh_t_latent = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report$pianka_o_cssif_sosh_t_latent,
                                                                       sosh_index_t_latent = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report$sosh_index_t_latent,
                                                                       pianka_o_cssif_comu_t_latent = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report$pianka_o_cssif_comu_t_latent,
                                                                       comu_index_t_latent = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report$comu_index_t_latent,
                                                                       pianka_o_cssif_prey_field_t_latent = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report$pianka_o_cssif_prey_field_t_latent,
                                                                       prey_field_index_t_latent = SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_report$prey_field_index_t_latent)

SRF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                         linear_predictor = 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_0","estimate"] +
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_sosh_t_latent +
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE[SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                           SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_comu_t_latent)

SRF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions$linear_predictor)

SRF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions$var <- SRF_06_2_seabird_prey_cssif_only_SAR_SEcov_pred_var

SRF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions

SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions %>% 
  left_join(SRF_subset_SAR, by = "run_year") -> SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp


SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp_plot <- ggplot(SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp_plot.png"), SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### UCSF - Yearlings x Seabirds ####
UCSF_subset %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> UCSF_subset_SAR

UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$variance <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$std_error^2

## variance of beta_0 parameter
var_beta_0 <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0", "variance"]

## csyif index of abundance
# expected value of beta_csyif index parameter
est_beta_csyif_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_csyif_index", "estimate"]

# variance of beta_csyif index parameter
var_beta_csyif_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_csyif_index", "variance"]

# expected value of csyif_index of abundance predictor
est_csyif_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$csyif_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$csyif_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$csyif_index_of_abundance)

# variance of csyif index of abundance predictor
var_csyif_index_of_abundance <- sqrt(diag(csyif_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$csyif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (UCSF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t - mean(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t))/sd(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_csyif_prey_field_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t))^2))


## sosh index
# expected value of beta_sosh index parameter
est_beta_sosh_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_sosh_index", "estimate"]

# variance of beta_sosh index parameter
var_beta_sosh_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_sosh_index", "variance"]

# expected value of sosh_index of abundance predictor
est_sosh_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance)

# variance of sosh index of abundance predictor
var_sosh_index_of_abundance <- sqrt(diag(sosh_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance))^2))

## sosh overlap
# expected value of beta_sosh parameter
est_beta_ov_sosh <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_sosh", "estimate"]

# variance of beta_sosh parameter
var_beta_ov_sosh <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_sosh", "variance"]

# expected value of sosh overlap predictor
est_sosh_overlap <- (UCSF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t - mean(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t))/sd(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t)

# variance of sosh index of abundance predictor
var_sosh_overlap <- sqrt(diag(pianka_o_csyif_sosh_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t))^2))

## comu index
# expected value of beta_comu index parameter
est_beta_comu_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_comu_index", "estimate"]

# variance of beta_comu index parameter
var_beta_comu_index <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_comu_index", "variance"]

# expected value of comu_index of abundance predictor
est_comu_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$comu_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$comu_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$comu_index_of_abundance)

# variance of comu index of abundance predictor
var_comu_index_of_abundance <- sqrt(diag(comu_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$comu_index_of_abundance))^2))

## comu overlap
# expected value of beta_comu parameter
est_beta_ov_comu <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_comu", "estimate"]

# variance of beta_comu parameter
var_beta_ov_comu <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_comu", "variance"]

# expected value of comu overlap predictor
est_comu_overlap <- (UCSF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t - mean(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t))/sd(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t)

# variance of comu index of abundance predictor
var_comu_overlap <- sqrt(diag(pianka_o_csyif_comu_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$pianka_o_csyif_comu_t))^2))










#### Calculate and plot total variance by year
UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
           est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                        csyif_index_t_latent = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report$csyif_index_t_latent,
                                                                        pianka_o_csyif_sosh_t_latent = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report$pianka_o_csyif_sosh_t_latent,
                                                                        sosh_index_t_latent = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report$sosh_index_t_latent,
                                                                        pianka_o_csyif_comu_t_latent = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report$pianka_o_csyif_comu_t_latent,
                                                                        comu_index_t_latent = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report$comu_index_t_latent,
                                                                        pianka_o_csyif_prey_field_t_latent = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report$pianka_o_csyif_prey_field_t_latent,
                                                                        prey_field_index_t_latent = UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_report$prey_field_index_t_latent)

UCSF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                          linear_predictor = 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0","estimate"] +
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_sosh_t_latent +
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_comu_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_comu_t_latent)

UCSF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions$linear_predictor)

# join the analytical solution with the MLE predictions

UCSF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions$var <- UCSF_06_2_seabird_prey_csyif_only_SAR_SEcov_pred_var

UCSF_06_2_seabird_prey_csyif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions

UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp


UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp_plot <- ggplot(UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp_plot.png"), UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### UCSF - Subyearlings x Seabirds ####
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$variance <-UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$std_error^2


## variance of beta_0 parameter
var_beta_0 <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0", "variance"]

## cssif index of abundance
# expected value of beta_cssif index parameter
est_beta_cssif_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_cssif_index", "estimate"]

# variance of beta_cssif index parameter
var_beta_cssif_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_cssif_index", "variance"]

# expected value of cssif_index of abundance predictor
est_cssif_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$cssif_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$cssif_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$cssif_index_of_abundance)

# variance of cssif index of abundance predictor
var_cssif_index_of_abundance <- sqrt(diag(cssif_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$cssif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (UCSF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t - mean(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t))/sd(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_cssif_prey_field_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t))^2))


## sosh index
# expected value of beta_sosh index parameter
est_beta_sosh_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_sosh_index", "estimate"]

# variance of beta_sosh index parameter
var_beta_sosh_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_sosh_index", "variance"]

# expected value of sosh_index of abundance predictor
est_sosh_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance)

# variance of sosh index of abundance predictor
var_sosh_index_of_abundance <- sqrt(diag(sosh_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$sosh_index_of_abundance))^2))

## sosh overlap
# expected value of beta_sosh parameter
est_beta_ov_sosh <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_sosh", "estimate"]

# variance of beta_sosh parameter
var_beta_ov_sosh <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_sosh", "variance"]

# expected value of sosh overlap predictor
est_sosh_overlap <- (UCSF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t - mean(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t))/sd(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t)

# variance of sosh index of abundance predictor
var_sosh_overlap <- sqrt(diag(pianka_o_cssif_sosh_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t))^2))

## comu index
# expected value of beta_comu index parameter
est_beta_comu_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_comu_index", "estimate"]

# variance of beta_comu index parameter
var_beta_comu_index <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_comu_index", "variance"]

# expected value of comu_index of abundance predictor
est_comu_index_of_abundance <- (UCSF_v1_seabird_prey_predictors$comu_index_of_abundance - mean(UCSF_v1_seabird_prey_predictors$comu_index_of_abundance))/sd(UCSF_v1_seabird_prey_predictors$comu_index_of_abundance)

# variance of comu index of abundance predictor
var_comu_index_of_abundance <- sqrt(diag(comu_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$comu_index_of_abundance))^2))

## comu overlap
# expected value of beta_comu parameter
est_beta_ov_comu <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_comu", "estimate"]

# variance of beta_comu parameter
var_beta_ov_comu <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_comu", "variance"]

# expected value of comu overlap predictor
est_comu_overlap <- (UCSF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t - mean(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t))/sd(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t)

# variance of comu index of abundance predictor
var_comu_overlap <- sqrt(diag(pianka_o_cssif_comu_cov_matrix_common_years_seabird_prey * (1/sd(UCSF_v1_seabird_prey_predictors$pianka_o_cssif_comu_t))^2))

#### Calculate and plot total variance by year
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
           est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                        cssif_index_t_latent = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report$cssif_index_t_latent,
                                                                        pianka_o_cssif_sosh_t_latent = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report$pianka_o_cssif_sosh_t_latent,
                                                                        sosh_index_t_latent = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report$sosh_index_t_latent,
                                                                        pianka_o_cssif_comu_t_latent = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report$pianka_o_cssif_comu_t_latent,
                                                                        comu_index_t_latent = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report$comu_index_t_latent,
                                                                        pianka_o_cssif_prey_field_t_latent = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report$pianka_o_cssif_prey_field_t_latent,
                                                                        prey_field_index_t_latent = UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_report$prey_field_index_t_latent)

UCSF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                          linear_predictor = 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0","estimate"] +
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_sosh_t_latent +
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_comu_index","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                            UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_comu_t_latent)

UCSF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions$linear_predictor)

UCSF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions$var <- UCSF_06_2_seabird_prey_cssif_only_SAR_SEcov_pred_var

UCSF_06_2_seabird_prey_cssif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions

UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp


UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp_plot <- ggplot(UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp_plot.png"), UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

#### SRF - Yearlings x Hake ####
SRF_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_subset_SAR

SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$variance <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$std_error^2

## variance of beta_0 parameter
var_beta_0 <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_0", "variance"]

## csyif index of abundance
# expected value of beta_csyif index parameter
est_beta_csyif_index <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index", "estimate"]

# variance of beta_csyif index parameter
var_beta_csyif_index <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index", "variance"]

# expected value of csyif_index of abundance predictor
est_csyif_index_of_abundance <- (SRF_v1_hake_prey_predictors$csyif_index_of_abundance - mean(SRF_v1_hake_prey_predictors$csyif_index_of_abundance))/sd(SRF_v1_hake_prey_predictors$csyif_index_of_abundance)

# variance of csyif index of abundance predictor
var_csyif_index_of_abundance <- sqrt(diag(csyif_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$csyif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (SRF_v1_hake_prey_predictors$prey_field_index_of_abundance - mean(SRF_v1_hake_prey_predictors$prey_field_index_of_abundance))/sd(SRF_v1_hake_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (SRF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t - mean(SRF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t))/sd(SRF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_csyif_prey_field_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t))^2))


## hake index
# expected value of beta_hake index parameter
est_beta_hake_index <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_hake_index", "estimate"]

# variance of beta_hake index parameter
var_beta_hake_index <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_hake_index", "variance"]

# expected value of hake_index of abundance predictor
est_hake_index_of_abundance <- (SRF_v1_hake_prey_predictors$hake_index_of_abundance - mean(SRF_v1_hake_prey_predictors$hake_index_of_abundance))/sd(SRF_v1_hake_prey_predictors$hake_index_of_abundance)

# variance of hake index of abundance predictor
var_hake_index_of_abundance <- sqrt(diag(hake_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$hake_index_of_abundance))^2))

## hake overlap
# expected value of beta_hake parameter
est_beta_ov_hake <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_hake", "estimate"]

# variance of beta_hake parameter
var_beta_ov_hake <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_hake", "variance"]

# expected value of hake overlap predictor
est_hake_overlap <- (SRF_v1_hake_prey_predictors$pianka_o_csyif_hake_t - mean(SRF_v1_hake_prey_predictors$pianka_o_csyif_hake_t))/sd(SRF_v1_hake_prey_predictors$pianka_o_csyif_hake_t)

# variance of hake index of abundance predictor
var_hake_overlap <- sqrt(diag(pianka_o_csyif_hake_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$pianka_o_csyif_hake_t))^2))


#### Calculate and plot total variance by year
SRF_06_3_hake_prey_csyif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
           est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
SRF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_hake_prey,
                                                                    csyif_index_t_latent = SRF_06_3_hake_prey_csyif_only_SAR_SEcov_report$csyif_index_t_latent,
                                                                    pianka_o_csyif_hake_t_latent = SRF_06_3_hake_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_hake_t_latent,
                                                                    hake_index_t_latent = SRF_06_3_hake_prey_csyif_only_SAR_SEcov_report$hake_index_t_latent,
                                                                    pianka_o_csyif_prey_field_t_latent = SRF_06_3_hake_prey_csyif_only_SAR_SEcov_report$pianka_o_csyif_prey_field_t_latent,
                                                                    prey_field_index_t_latent = SRF_06_3_hake_prey_csyif_only_SAR_SEcov_report$prey_field_index_t_latent)

SRF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                      linear_predictor = 
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_0","estimate"] +
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_csyif_only_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                        SRF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_hake_t_latent)

SRF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions$linear_predictor)

SRF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions$var <- SRF_06_3_hake_prey_csyif_only_SAR_SEcov_pred_var

SRF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_06_3_hake_prey_csyif_only_SAR_model_predictions

SRF_06_3_hake_prey_csyif_only_SAR_model_predictions %>% 
  left_join(SRF_subset_SAR, by = "run_year") -> SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp


SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp_plot <- ggplot(SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp_plot.png"), SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### SRF - Subyearlings x Hake ####
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$variance <-SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$std_error^2

## variance of beta_0 parameter
var_beta_0 <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_0", "variance"]

## cssif index of abundance
# expected value of beta_cssif index parameter
est_beta_cssif_index <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_cssif_index", "estimate"]

# variance of beta_cssif index parameter
var_beta_cssif_index <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_cssif_index", "variance"]

# expected value of cssif_index of abundance predictor
est_cssif_index_of_abundance <- (SRF_v1_hake_prey_predictors$cssif_index_of_abundance - mean(SRF_v1_hake_prey_predictors$cssif_index_of_abundance))/sd(SRF_v1_hake_prey_predictors$cssif_index_of_abundance)

# variance of cssif index of abundance predictor
var_cssif_index_of_abundance <- sqrt(diag(cssif_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$cssif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (SRF_v1_hake_prey_predictors$prey_field_index_of_abundance - mean(SRF_v1_hake_prey_predictors$prey_field_index_of_abundance))/sd(SRF_v1_hake_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (SRF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t - mean(SRF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t))/sd(SRF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_cssif_prey_field_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t))^2))


## hake index
# expected value of beta_hake index parameter
est_beta_hake_index <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_hake_index", "estimate"]

# variance of beta_hake index parameter
var_beta_hake_index <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_hake_index", "variance"]

# expected value of hake_index of abundance predictor
est_hake_index_of_abundance <- (SRF_v1_hake_prey_predictors$hake_index_of_abundance - mean(SRF_v1_hake_prey_predictors$hake_index_of_abundance))/sd(SRF_v1_hake_prey_predictors$hake_index_of_abundance)

# variance of hake index of abundance predictor
var_hake_index_of_abundance <- sqrt(diag(hake_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$hake_index_of_abundance))^2))

## hake overlap
# expected value of beta_hake parameter
est_beta_ov_hake <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_hake", "estimate"]

# variance of beta_hake parameter
var_beta_ov_hake <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_hake", "variance"]

# expected value of hake overlap predictor
est_hake_overlap <- (SRF_v1_hake_prey_predictors$pianka_o_cssif_hake_t - mean(SRF_v1_hake_prey_predictors$pianka_o_cssif_hake_t))/sd(SRF_v1_hake_prey_predictors$pianka_o_cssif_hake_t)

# variance of hake index of abundance predictor
var_hake_overlap <- sqrt(diag(pianka_o_cssif_hake_cov_matrix_common_years_hake_prey * (1/sd(SRF_v1_hake_prey_predictors$pianka_o_cssif_hake_t))^2))


#### Calculate and plot total variance by year
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
           est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
SRF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_hake_prey,
                                                                    cssif_index_t_latent = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_report$cssif_index_t_latent,
                                                                    pianka_o_cssif_hake_t_latent = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_report$pianka_o_cssif_hake_t_latent,
                                                                    hake_index_t_latent = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_report$hake_index_t_latent,
                                                                    pianka_o_cssif_prey_field_t_latent = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_report$pianka_o_cssif_prey_field_t_latent,
                                                                    prey_field_index_t_latent = SRF_06_3_hake_prey_cssif_only_SAR_SEcov_report$prey_field_index_t_latent)

SRF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                      linear_predictor = 
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_0","estimate"] +
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE[SRF_06_3_hake_prey_cssif_only_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                        SRF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_hake_t_latent)

SRF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions$linear_predictor)

SRF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions$var <- SRF_06_3_hake_prey_cssif_only_SAR_SEcov_pred_var

SRF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_06_3_hake_prey_cssif_only_SAR_model_predictions

SRF_06_3_hake_prey_cssif_only_SAR_model_predictions %>% 
  left_join(SRF_subset_SAR, by = "run_year") -> SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp


SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp_plot <- ggplot(SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp_plot.png"), SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### UCSF - Yearlings x Hake ####
UCSF_subset %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> UCSF_subset_SAR

UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$variance <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$std_error^2

## variance of beta_0 parameter
var_beta_0 <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0", "variance"]

## csyif index of abundance
# expected value of beta_csyif index parameter
est_beta_csyif_index <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_csyif_index", "estimate"]

# variance of beta_csyif index parameter
var_beta_csyif_index <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_csyif_index", "variance"]

# expected value of csyif_index of abundance predictor
est_csyif_index_of_abundance <- (UCSF_v1_hake_prey_predictors$csyif_index_of_abundance - mean(UCSF_v1_hake_prey_predictors$csyif_index_of_abundance))/sd(UCSF_v1_hake_prey_predictors$csyif_index_of_abundance)

# variance of csyif index of abundance predictor
var_csyif_index_of_abundance <- sqrt(diag(csyif_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$csyif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance - mean(UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance))/sd(UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (UCSF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t - mean(UCSF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t))/sd(UCSF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_csyif_prey_field_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t))^2))


## hake index
# expected value of beta_hake index parameter
est_beta_hake_index <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_hake_index", "estimate"]

# variance of beta_hake index parameter
var_beta_hake_index <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_hake_index", "variance"]

# expected value of hake_index of abundance predictor
est_hake_index_of_abundance <- (UCSF_v1_hake_prey_predictors$hake_index_of_abundance - mean(UCSF_v1_hake_prey_predictors$hake_index_of_abundance))/sd(UCSF_v1_hake_prey_predictors$hake_index_of_abundance)

# variance of hake index of abundance predictor
var_hake_index_of_abundance <- sqrt(diag(hake_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$hake_index_of_abundance))^2))

## hake overlap
# expected value of beta_hake parameter
est_beta_ov_hake <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_hake", "estimate"]

# variance of beta_hake parameter
var_beta_ov_hake <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_hake", "variance"]

# expected value of hake overlap predictor
est_hake_overlap <- (UCSF_v1_hake_prey_predictors$pianka_o_csyif_hake_t - mean(UCSF_v1_hake_prey_predictors$pianka_o_csyif_hake_t))/sd(UCSF_v1_hake_prey_predictors$pianka_o_csyif_hake_t)

# variance of hake index of abundance predictor
var_hake_overlap <- sqrt(diag(pianka_o_csyif_hake_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$pianka_o_csyif_hake_t))^2))


#### Calculate and plot total variance by year
UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
           est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_hake_prey,
                                                                     csyif_index_t_latent = UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_report$csyif_index_t_latent,
                                                                     pianka_o_csyif_hake_t_latent = UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_report$pianka_o_csyif_hake_t_latent,
                                                                     hake_index_t_latent = UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_report$hake_index_t_latent,
                                                                     pianka_o_csyif_prey_field_t_latent = UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_report$pianka_o_csyif_prey_field_t_latent,
                                                                     prey_field_index_t_latent = UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_report$prey_field_index_t_latent)

UCSF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                       linear_predictor = 
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0","estimate"] +
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_hake_index","estimate"] * 
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                         UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_predictors_df$pianka_o_csyif_hake_t_latent)

UCSF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions$linear_predictor)

# join the analytical solution with the MLE predictions

UCSF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions$var <- UCSF_06_3_hake_prey_csyif_only_SAR_SEcov_pred_var

UCSF_06_3_hake_prey_csyif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions

UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp


UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp_plot <- ggplot(UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp_plot.png"), UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### UCSF - Subyearlings x Hake ####
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$variance <-UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$std_error^2


## variance of beta_0 parameter
var_beta_0 <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0", "variance"]

## cssif index of abundance
# expected value of beta_cssif index parameter
est_beta_cssif_index <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_cssif_index", "estimate"]

# variance of beta_cssif index parameter
var_beta_cssif_index <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_cssif_index", "variance"]

# expected value of cssif_index of abundance predictor
est_cssif_index_of_abundance <- (UCSF_v1_hake_prey_predictors$cssif_index_of_abundance - mean(UCSF_v1_hake_prey_predictors$cssif_index_of_abundance))/sd(UCSF_v1_hake_prey_predictors$cssif_index_of_abundance)

# variance of cssif index of abundance predictor
var_cssif_index_of_abundance <- sqrt(diag(cssif_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$cssif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance - mean(UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance))/sd(UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (UCSF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t - mean(UCSF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t))/sd(UCSF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_cssif_prey_field_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t))^2))


## hake index
# expected value of beta_hake index parameter
est_beta_hake_index <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_hake_index", "estimate"]

# variance of beta_hake index parameter
var_beta_hake_index <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_hake_index", "variance"]

# expected value of hake_index of abundance predictor
est_hake_index_of_abundance <- (UCSF_v1_hake_prey_predictors$hake_index_of_abundance - mean(UCSF_v1_hake_prey_predictors$hake_index_of_abundance))/sd(UCSF_v1_hake_prey_predictors$hake_index_of_abundance)

# variance of hake index of abundance predictor
var_hake_index_of_abundance <- sqrt(diag(hake_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$hake_index_of_abundance))^2))

## hake overlap
# expected value of beta_hake parameter
est_beta_ov_hake <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_hake", "estimate"]

# variance of beta_hake parameter
var_beta_ov_hake <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_hake", "variance"]

# expected value of hake overlap predictor
est_hake_overlap <- (UCSF_v1_hake_prey_predictors$pianka_o_cssif_hake_t - mean(UCSF_v1_hake_prey_predictors$pianka_o_cssif_hake_t))/sd(UCSF_v1_hake_prey_predictors$pianka_o_cssif_hake_t)

# variance of hake index of abundance predictor
var_hake_overlap <- sqrt(diag(pianka_o_cssif_hake_cov_matrix_common_years_hake_prey * (1/sd(UCSF_v1_hake_prey_predictors$pianka_o_cssif_hake_t))^2))

#### Calculate and plot total variance by year
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_pred_var <- var_beta_0 + 
  var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
           est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
# get predictors in a df
UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df <- data.frame(year = common_years_hake_prey,
                                                                     cssif_index_t_latent = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_report$cssif_index_t_latent,
                                                                     pianka_o_cssif_hake_t_latent = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_report$pianka_o_cssif_hake_t_latent,
                                                                     hake_index_t_latent = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_report$hake_index_t_latent,
                                                                     pianka_o_cssif_prey_field_t_latent = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_report$pianka_o_cssif_prey_field_t_latent,
                                                                     prey_field_index_t_latent = UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_report$prey_field_index_t_latent)

UCSF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                       linear_predictor = 
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_0","estimate"] +
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_hake_index","estimate"] * 
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE[UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_no_transport_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                         UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_predictors_df$pianka_o_cssif_hake_t_latent)

UCSF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions$linear_predictor)

UCSF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions$var <- UCSF_06_3_hake_prey_cssif_only_SAR_SEcov_pred_var

UCSF_06_3_hake_prey_cssif_only_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions

UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp


UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp_plot <- ggplot(UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical")) +
  geom_point(data = UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted")) +
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

ggsave(here::here("figures", "paper_figures", "SAR_prediction_plots", "UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp_plot.png"), UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### Combine figure ####

SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp$model <- "Snake River Fall - Yearlings x Seabirds"
SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp$model <- "Snake River Fall - Yearlings x Hake"
SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp$model  <- "Snake River Fall - Subyearlings x Seabirds"
SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp$model  <- "Snake River Fall - Subyearlings x Hake"

UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp$model <- "Upper Columbia Summer/Fall - Yearlings x Seabirds"
UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp$model <- "Upper Columbia Summer/Fall - Yearlings x Hake"
UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp$model <- "Upper Columbia Summer/Fall - Subyearlings x Seabirds"
UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp$model <- "Upper Columbia Summer/Fall - Subyearlings x Hake"

SRF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp %>% 
  bind_rows(SRF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp) %>% 
  bind_rows(UCSF_06_2_seabird_prey_csyif_only_SAR_model_predictions_comp) %>% 
  bind_rows(UCSF_06_2_seabird_prey_cssif_only_SAR_model_predictions_comp)  -> SAR_seabird_model_predictions_comp
  
SRF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp %>%
  bind_rows(SRF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp) %>%
  bind_rows(UCSF_06_3_hake_prey_csyif_only_SAR_model_predictions_comp) %>%
  bind_rows(UCSF_06_3_hake_prey_cssif_only_SAR_model_predictions_comp) -> SAR_hake_model_predictions_comp

SAR_seabird_model_predictions_comp %>% 
  bind_rows(SAR_hake_model_predictions_comp) -> SAR_model_predictions_comp

SAR_model_predictions_comp$model <- factor(SAR_model_predictions_comp$model, 
                                           levels = c("Snake River Fall - Subyearlings x Seabirds",
                                                      "Snake River Fall - Subyearlings x Hake",
                                                      "Snake River Fall - Yearlings x Seabirds",
                                                      "Snake River Fall - Yearlings x Hake",
                                                      "Upper Columbia Summer/Fall - Subyearlings x Seabirds",
                                                      "Upper Columbia Summer/Fall - Subyearlings x Hake",
                                                      "Upper Columbia Summer/Fall - Yearlings x Seabirds",
                                                      "Upper Columbia Summer/Fall - Yearlings x Hake"))
labels_df <- data.frame(model = c("Snake River Fall - Subyearlings x Seabirds",
                                  "Snake River Fall - Subyearlings x Hake",
                                  "Snake River Fall - Yearlings x Seabirds",
                                  "Snake River Fall - Yearlings x Hake",
                                  "Upper Columbia Summer/Fall - Subyearlings x Seabirds",
                                  "Upper Columbia Summer/Fall - Subyearlings x Hake",
                                  "Upper Columbia Summer/Fall - Yearlings x Seabirds",
                                  "Upper Columbia Summer/Fall - Yearlings x Hake"),
                        label = c("(A) Snake River Fall - Subyearlings x Seabirds",
                                  "(B) Snake River Fall - Subyearlings x Hake",
                                  "(C) Snake River Fall - Yearlings x Seabirds",
                                  "(D) Snake River Fall - Yearlings x Hake",
                                  "(E) Upper Columbia Summer/Fall - Subyearlings x Seabirds",
                                  "(F) Upper Columbia Summer/Fall - Subyearlings x Hake",
                                  "(G) Upper Columbia Summer/Fall - Yearlings x Seabirds",
                                  "(H) Upper Columbia Summer/Fall - Yearlings x Hake"))

SAR_model_predictions_comp %>% 
  left_join(labels_df, by = "model") -> SAR_model_predictions_comp

fig5_SAR_model_predictions_comp_plot <- ggplot(SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  geom_errorbar(data = SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1) +
  geom_point(aes(shape = "Empirical"), size = 2) +
  geom_point(data = SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted"), size = 3) +
  scale_shape_manual(name = NULL, values = c("Empirical" = 19, "Predicted" = 8)) +
  xlab("Run Year") +
  ylab("Marine Survival") +
  # coord_cartesian(ylim=c(0, 0.1)) +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.position = "bottom",
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, hjust = 0),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))+
    facet_wrap(~label, ncol = 2)


ggsave(here::here("figures", "paper_figures", "fig5_SAR_model_predictions_comp_plot.png"), 
       fig5_SAR_model_predictions_comp_plot,  
       height = 12, width = 12)
