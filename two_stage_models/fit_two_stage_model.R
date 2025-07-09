## Two stage model
# last updated: 2025-07-08

# load libraries
library(tidyverse)
library(here)
library(TMB)

SAR_data <- read.csv(here::here("model_inputs", "chinook_det_hist.csv"))

# subset the SAR data to drop 1998 (we don't have bongo data in this year) and
# also drop 2022-2025 (since we don't have full adult returns for those years yet)
SAR_data <- subset(SAR_data, !(run_year %in% c(1998, 2022:2025)))

# source the scripts to prepare the data for the SDMs
source(here::here("SDM_models", "jsoes_bongo_SDM_prep.R"))
source(here::here("SDM_models", "jsoes_trawl_SDM_prep.R"))

#### prep model data ####


# extract SDM data
csyif <- subset(jsoes_long, species == "chinook_salmon_yearling_interior_fa")

# extract Snake River Fall Chinook from SAR data
SRF <- subset(SAR_data, run_name == "Fall" & group == "Snake River")

# generate the projection vector/matrix for covariates - this relies
# on the survey prediction grid from the earlier scripts
distance_projection_vector <- as.numeric(filter(survey_predict_grid, year == 2000)$dist_shore_scaled)

# for SST - take prediction grid, add SST with locations as rows and years as column
survey_predict_grid %>% 
  dplyr::select(-c(SST, dist_shore, dist_shore_scaled)) %>% 
  pivot_wider(names_from = year, values_from = SST_scaled) %>% 
  dplyr::select(-c(X, Y)) %>% 
  as.matrix() -> temperature_projection_matrix

# testing - let's see if rescaling the data for shrimp larvae would help
jsoes_bongo_shrimp_larvae %>% 
  mutate(total_sum_of_density_number_liter = total_sum_of_density_number_m3/1000) -> jsoes_bongo_shrimp_larvae


#### Stage 1: Fit the two SDMs ####

compile(here::here("two_stage_models", "SDMs_overlap.cpp"))
dyn.load(dynlib(here::here("two_stage_models", "SDMs_overlap")))

Data = list(
  ## Shared data between species A and species B
  "dist_g" = distance_projection_vector, "temp_gt" = temperature_projection_matrix, 
  "n_t" = length(unique(jsoes_bongo_shrimp_larvae$year)),
  
  ## Species A data (shrimp larvae)
  "D_j_a"= jsoes_bongo_shrimp_larvae$total_sum_of_density_number_liter, "t_j" = jsoes_bongo_shrimp_larvae$year - min(jsoes_bongo_shrimp_larvae$year),
  "temp_j" = jsoes_bongo_shrimp_larvae$SST_scaled, "dist_j" = jsoes_bongo_shrimp_larvae$dist_shore_scaled,
  "A_is"=bongo_A_is, "A_gs"=bongo_A_gs, "M0_a"=bongo_spde$c0, 
  "M1_a"=bongo_spde$g1, "M2_a"=bongo_spde$g2,
  
  ## Species B data (Yearling Interior Fall Chinook)
  "D_k_b"= csyif$n_per_km, "t_k" = csyif$year - min(csyif$year),
  "temp_k" = csyif$SST_scaled, "dist_k" = csyif$dist_shore_scaled,
  "B_is"=chinook_is, "B_gs"=chinook_gs, "M0_b"=chinook_spde$c0, "M1_b"=chinook_spde$g1, "M2_b"=chinook_spde$g2
)

Params <- list(
  ## SAR parameters
  "epsilon_t" = rep(0, length(unique(SRF$run_year))),
  "ln_sigma_t" = 0,
  
  ## Species A SDM parameters
  "beta_t_a"=rep(0, length(unique(jsoes_bongo_shrimp_larvae$year))),
  "beta_temp_a" = 0, "beta_dist_a" = 0,
  "ln_tau_omega_a"=0, "ln_tau_epsilon_a"=0,
  "ln_kappa_a"=0, "ln_phi_a" = 0,
  "finv_power_a" = 0,
  
  ## Species B SDM parameters
  "beta_t_b"=rep(0, length(unique(csyif$year))),
  "beta_temp_b" = 0, "beta_dist_b" = 0,
  "ln_tau_omega_b"=0, "ln_tau_epsilon_b"=0,
  "ln_kappa_b"=0, "ln_phi_b" = 0,
  "finv_power_b" = 0,
  
  # Species A SDM random effects
  "omega_s_a"=rnorm(nrow(bongo_spde$c0)),
  "epsilon_st_a"=matrix(0, nrow=nrow(bongo_spde$c0), ncol=Data$n_t),
  
  # Species B SDM random effects
  "omega_s_b"=rnorm(nrow(chinook_spde$c0)),
  "epsilon_st_b"=matrix(0, nrow=nrow(chinook_spde$c0), ncol=Data$n_t)
  
)

SDMs_overlap_Obj = MakeADFun( data=Data, 
                                    parameters=Params, 
                                    random= c("omega_s_a", "epsilon_st_a", "omega_s_b", "epsilon_st_b"),
                                    DLL = "SDMs_overlap")

# Optimize
SDMs_overlap_Opt = nlminb( start=SDMs_overlap_Obj$par, obj=SDMs_overlap_Obj$fn, grad=SDMs_overlap_Obj$gr )
# add getJointPrecision for index
SDMs_overlap_Opt$SD = sdreport( SDMs_overlap_Obj, bias.correct=TRUE, getJointPrecision = TRUE )
SDMs_overlap_report = SDMs_overlap_Obj$report()

# Extract the covariance for the index

## Use Jim's code
# Function to generate samples
sample_var <-
  function( obj, # TMB object
            var_name, # name of REPORTed variable
            mu, # Estimated mean
            prec, # Estimated precision
            n_samples = 500, # Number of samples
            fun1 = function(x,...) x, # optional samples transform
            fun2 = sd, # calculate fun2 across samples
            ... ){
    
    require(abind)
    # Take samples using precision matrix
    u_zr = rmvnorm_prec( mu=mu, prec=prec, n.sims=n_samples )
    # Calculate REPORTed variable for each sample
    for( rI in 1:n_samples ){
      Var = obj$report( par=u_zr[,rI] )[[var_name]]
      # Transform REPORTed variable using fun1
      Var = fun1( Var, ... )
      # Bind samples together into array
      if(is.vector(Var)) Var = as.array(Var)
      if(rI==1) Var_zr = Var
      if(rI>=2){
        Var_zr = abind( Var_zr, Var, along=length(dim(Var))+1 )
      }
    }
    # summarize across samples using fun2
    out = apply(Var_zr, MARGIN=1:(length(dim(Var_zr))-1), FUN=fun2)
    # Return value
    return( out )
  }

rmvnorm_prec <-
  function( mu, # estimated fixed and random effects
            prec, # estimated joint precision
            n.sims) {
    
    require(Matrix)
    # Simulate values
    z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    # Q = t(P) * L * t(L) * P
    L = Cholesky(prec, super=TRUE)
    # Calcualte t(P) * solve(t(L)) * z0 in two steps
    z = solve(L, z0, system = "Lt") # z = Lt^-1 * z
    z = solve(L, z, system = "Pt") # z = Pt    * z
    return(mu + as.matrix(z))
  }

SE_ov_ab_t = sample_var( obj=SDMs_overlap_Obj, var_name="ov_ab_t", mu=SDMs_overlap_Obj$env$last.par.best, prec=SDMs_overlap_Opt$SD$jointPrecision )





#### Stage 2: Fit the SAR model, using error-in-variables

compile(here::here("two_stage_models", "SAR_error_in_variables.cpp"))
dyn.load(dynlib(here::here("two_stage_models", "SAR_error_in_variables")))

Data = list(
  ## SAR data
  "z_i" = SRF$adult_det, "t_i" = SRF$run_year-min(SRF$run_year),
  
  ## data from SDM-derived overlap
  "SE_ov_ab_t" = SE_ov_ab_t, ov_ab_t = SDMs_overlap_report$ov_ab_t

)

Params <- list(
  ## SAR parameters
  "epsilon_t" = rep(0, length(unique(SRF$run_year))),
  "ln_sigma_t" = 0,
  "beta_ov_ab" = 0,
  
  ## error-in-variables parameters
  "ln_sigma_error" # this is the sigma for the unexplained residual error
  
)

SAR_error_in_variables_Obj = MakeADFun( data=Data, 
                              parameters=Params, 
                              random= c("epsilon_t"),
                              DLL = "SAR_error_in_variables")

# Optimize
SAR_error_in_variables_Opt = nlminb( start=SAR_error_in_variables_Obj$par, obj=SAR_error_in_variables_Obj$fn, grad=SAR_error_in_variables_Obj$gr )
SAR_error_in_variables_Opt$SD = sdreport( SAR_error_in_variables_Obj, bias.correct=TRUE )
SAR_error_in_variables_report = SAR_error_in_variables_Obj$report()