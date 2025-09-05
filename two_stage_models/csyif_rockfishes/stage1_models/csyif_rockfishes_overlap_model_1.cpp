// Using the best-fit SDMs for both CSYIF and rockfishes, estimate the overlap
// between the two within the survey area


#include <TMB.hpp>
#include <algorithm>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{

using namespace density;

// This model uses a lot of indexing to keep our taxa straight
// the subscript csyif denotes Yearling Interior Fall Chinook
// the subscript rf denotes Sebastes spp.
// the subscript i denotes samples from the JSOES trawl (for CSYIF)
// the subscript i=j denotes samples from the PWCC or PRS surveys (for rockfishes)

// Objective function
Type jnll = 0;

// #### SHARED DATA ####

DATA_INTEGER(n_t); // number of years in the dataset
// DATA_MATRIX(temp_gt); // temperature at each location in each year
// DATA_VECTOR(dist_g); // distance from shore at each location


// #### CSYIF SDM ####

// Data for csyif
DATA_VECTOR(D_i_csyif);  // density of csyif in measurement i
DATA_IVECTOR(t_i_csyif); // index for the year of measurement i (csyif)
DATA_VECTOR(weights_i_csyif); // optional weights for csyif - used for cAIC


// Projection matrices for csyif
DATA_SPARSE_MATRIX(A_is_csyif); // the projection matrix from jsoes trawl spde vertices to the jsoes trawl samples
DATA_SPARSE_MATRIX(A_gs_csyif); // the projection matrix from jsoes trawl spde vertices to the projection grid for the survey domain


// Parameters for csyif
PARAMETER_VECTOR( beta_t_csyif ); // fixed effects for year on density of csyif
PARAMETER( ln_tau_omega_csyif ); // Tau parameter for spatial effects (omega) for csyif
PARAMETER( ln_tau_epsilon_csyif ); // Tau parameter for spatiotemporal effects (epsilon) for csyif
PARAMETER( ln_kappa_csyif ); // Kappa term for SPDE for csyif
PARAMETER( ln_phi_csyif ); // phi term in tweedie for csyif
PARAMETER( finv_power_csyif ); // power parameter in tweedie for csyif

PARAMETER_VECTOR( ln_H_input_csyif ); // anisotropy input for csyif

// Random effects for csyif
PARAMETER_VECTOR( omega_s_csyif ); // vector of spatial random effects for csyif
PARAMETER_MATRIX( epsilon_st_csyif ); // vector of spatiotemporal random effects for csyif

int n_i_csyif = A_is_csyif.rows(); // number of csyif samples
int n_g_csyif = A_gs_csyif.rows(); // number of units in the projection grid for csyif

// priors for anisotropy parameters for csyif
Type ln_H_0_csyif_mean = 0.0;
Type ln_H_0_csyif_sd = 0.5;
Type ln_H_1_csyif_mean = 0.0;
Type ln_H_1_csyif_sd = 0.5;

jnll -= dnorm(ln_H_input_csyif(0), ln_H_0_csyif_mean, ln_H_0_csyif_sd, true); // Northings anisotropy
jnll -= dnorm(ln_H_input_csyif(1), ln_H_1_csyif_mean, ln_H_1_csyif_sd, true); // Anisotropic correlation

// Anisotropy elements for csyif
matrix<Type> H_csyif( 2, 2 );
H_csyif(0,0) = exp(ln_H_input_csyif(0));
H_csyif(1,0) = ln_H_input_csyif(1);
H_csyif(0,1) = ln_H_input_csyif(1);
H_csyif(1,1) = (1+ln_H_input_csyif(1)*ln_H_input_csyif(1)) / exp(ln_H_input_csyif(0));

// implement anisotropy for csyif
Eigen::SparseMatrix<Type> Q_csyif;
// Using INLA
DATA_STRUCT( spatial_list_csyif, R_inla::spde_aniso_t );
// Build precision matrix for csyif
Q_csyif = R_inla::Q_spde( spatial_list_csyif, exp(ln_kappa_csyif), H_csyif );

// Derived quantities
// SDM derived quantities for csyif
Type Range_csyif = sqrt(8) / exp( ln_kappa_csyif );
Type SigmaO_csyif = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_csyif) * exp(2*ln_kappa_csyif));
Type SigmaE_csyif = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_csyif) * exp(2*ln_kappa_csyif));

// Probability of random effects for csyif
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_csyif), 1/exp(ln_tau_omega_csyif) )( omega_s_csyif );
// spatio-temporal random effect - scaled by ln_tau_epsilon
for( int t=0; t<n_t; t++){
  jnll += SCALE( GMRF(Q_csyif), 1/exp(ln_tau_epsilon_csyif) )( epsilon_st_csyif.col(t) );
}

// SDM projections for csyif
vector<Type> omega_g_csyif( n_g_csyif );
omega_g_csyif = A_gs_csyif * omega_s_csyif;
matrix<Type> epsilon_gt_csyif( n_g_csyif, n_t );
epsilon_gt_csyif = A_gs_csyif * epsilon_st_csyif;

// Probability of data conditional on random effects for csyif
vector<Type> omega_i_csyif( n_i_csyif );
omega_i_csyif = A_is_csyif * omega_s_csyif;
matrix<Type> epsilon_it_csyif( n_i_csyif, n_t );
epsilon_it_csyif = A_is_csyif * epsilon_st_csyif;

vector<Type> dhat_i_csyif( n_i_csyif );
dhat_i_csyif.setZero();


for( int i=0; i<D_i_csyif.size(); i++){
  dhat_i_csyif(i) = exp( beta_t_csyif(t_i_csyif(i)) + omega_i_csyif(i) + epsilon_it_csyif(i, t_i_csyif(i)));
  jnll -= dtweedie( D_i_csyif(i), dhat_i_csyif(i), exp(ln_phi_csyif), Type(1.0)+invlogit(finv_power_csyif), true ) * weights_i_csyif(i);
}



// SDM: modeled density in the projection grid for csyif
// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_csyif( n_g_csyif, n_t );
for( int t=0; t<n_t; t++){
  for( int g=0; g<n_g_csyif; g++){
    ln_d_gt_csyif(g,t) = beta_t_csyif(t) + omega_g_csyif(g) + epsilon_gt_csyif(g,t);
  }
}

// Reporting: csyif
REPORT( dhat_i_csyif );
REPORT( beta_t_csyif );
REPORT( ln_tau_omega_csyif );
REPORT( ln_tau_epsilon_csyif );
REPORT( ln_kappa_csyif );
REPORT( ln_phi_csyif );
REPORT( finv_power_csyif );
REPORT(SigmaO_csyif);
REPORT(SigmaE_csyif);
REPORT(Range_csyif);
REPORT(ln_d_gt_csyif);
REPORT( H_csyif );



// #### ROCKFISH SDM ####

// Data for rf
DATA_VECTOR(D_i_rf);  // density of rf in measurement i
DATA_IVECTOR(t_i_rf); // index for the year of measurement i (rf)
DATA_VECTOR(weights_i_rf); // optional weights for rf - used for cAIC


// Projection matrices for rf
DATA_SPARSE_MATRIX(A_is_rf); // the projection matrix from jsoes trawl spde vertices to the jsoes trawl samples
DATA_SPARSE_MATRIX(A_gs_rf); // the projection matrix from jsoes trawl spde vertices to the projection grid for the survey domain


// Parameters for rf
PARAMETER_VECTOR( beta_t_rf ); // fixed effects for year on density of rf
PARAMETER( ln_tau_omega_rf ); // Tau parameter for spatial effects (omega) for rf
PARAMETER( ln_tau_epsilon_rf ); // Tau parameter for spatiotemporal effects (epsilon) for rf
PARAMETER( ln_kappa_rf ); // Kappa term for SPDE for rf
PARAMETER( ln_phi_rf ); // phi term in tweedie for rf
PARAMETER( finv_power_rf ); // power parameter in tweedie for rf

PARAMETER_VECTOR( ln_H_input_rf ); // anisotropy input for rf

// Random effects for rf
PARAMETER_VECTOR( omega_s_rf ); // vector of spatial random effects for rf
PARAMETER_MATRIX( epsilon_st_rf ); // vector of spatiotemporal random effects for rf

int n_i_rf = A_is_rf.rows(); // number of rf samples
int n_g_rf = A_gs_rf.rows(); // number of units in the projection grid for rf

// priors for anisotropy parameters for rf
Type ln_H_0_rf_mean = 0.0;
Type ln_H_0_rf_sd = 0.5;
Type ln_H_1_rf_mean = 0.0;
Type ln_H_1_rf_sd = 0.5;

jnll -= dnorm(ln_H_input_rf(0), ln_H_0_rf_mean, ln_H_0_rf_sd, true); // Northings anisotropy
jnll -= dnorm(ln_H_input_rf(1), ln_H_1_rf_mean, ln_H_1_rf_sd, true); // Anisotropic correlation

// Anisotropy elements for rf
matrix<Type> H_rf( 2, 2 );
H_rf(0,0) = exp(ln_H_input_rf(0));
H_rf(1,0) = ln_H_input_rf(1);
H_rf(0,1) = ln_H_input_rf(1);
H_rf(1,1) = (1+ln_H_input_rf(1)*ln_H_input_rf(1)) / exp(ln_H_input_rf(0));

// implement anisotropy for rf
Eigen::SparseMatrix<Type> Q_rf;
// Using INLA
DATA_STRUCT( spatial_list_rf, R_inla::spde_aniso_t );
// Build precision matrix for rf
Q_rf = R_inla::Q_spde( spatial_list_rf, exp(ln_kappa_rf), H_rf );

// Derived quantities
// SDM derived quantities for rf
Type Range_rf = sqrt(8) / exp( ln_kappa_rf );
Type SigmaO_rf = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_rf) * exp(2*ln_kappa_rf));
Type SigmaE_rf = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_rf) * exp(2*ln_kappa_rf));

// Probability of random effects for rf
// spatial random effect - scaled by ln_tau_omega
jnll += SCALE( GMRF(Q_rf), 1/exp(ln_tau_omega_rf) )( omega_s_rf );
// spatio-temporal random effect - scaled by ln_tau_epsilon
for( int t=0; t<n_t; t++){
  jnll += SCALE( GMRF(Q_rf), 1/exp(ln_tau_epsilon_rf) )( epsilon_st_rf.col(t) );
}

// SDM projections for rf
vector<Type> omega_g_rf( n_g_rf );
omega_g_rf = A_gs_rf * omega_s_rf;
matrix<Type> epsilon_gt_rf( n_g_rf, n_t );
epsilon_gt_rf = A_gs_rf * epsilon_st_rf;

// Probability of data conditional on random effects for rf
vector<Type> omega_i_rf( n_i_rf );
omega_i_rf = A_is_rf * omega_s_rf;
matrix<Type> epsilon_it_rf( n_i_rf, n_t );
epsilon_it_rf = A_is_rf * epsilon_st_rf;

vector<Type> dhat_i_rf( n_i_rf );
dhat_i_rf.setZero();


for( int i=0; i<D_i_rf.size(); i++){
  dhat_i_rf(i) = exp( beta_t_rf(t_i_rf(i)) + omega_i_rf(i) + epsilon_it_rf(i, t_i_rf(i)));
  jnll -= dtweedie( D_i_rf(i), dhat_i_rf(i), exp(ln_phi_rf), Type(1.0)+invlogit(finv_power_rf), true ) * weights_i_rf(i);
}



// SDM: modeled density in the projection grid for rf
// Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
array<Type> ln_d_gt_rf( n_g_rf, n_t );
for( int t=0; t<n_t; t++){
  for( int g=0; g<n_g_rf; g++){
    ln_d_gt_rf(g,t) = beta_t_rf(t) + omega_g_rf(g) + epsilon_gt_rf(g,t);
  }
}

// Reporting: rf
REPORT( dhat_i_rf );
REPORT( beta_t_rf );
REPORT( ln_tau_omega_rf );
REPORT( ln_tau_epsilon_rf );
REPORT( ln_kappa_rf );
REPORT( ln_phi_rf );
REPORT( finv_power_rf );
REPORT(SigmaO_rf);
REPORT(SigmaE_rf);
REPORT(Range_rf);
REPORT(ln_d_gt_rf);
REPORT( H_rf );


// #### CALCULATE OVERLAP ####

// Use the modeled density in the projection grid to calculate biomass-weighted overlap
// Use a nested for loop to loop through all years and all grid cells
// note that n_g_csyif and n_g_rf should have the same dimensions and therefore we can use either to construct the array

// construct two arrays: one to store the numerator and one to store the denominator for the index
array<Type> ov_gt_rf_csyif( n_g_csyif, n_t );
array<Type> scaled_dens_gt_csyif( n_g_csyif, n_t );

// construct a matrix to store the overlap index in each year
vector<Type> ov_rf_csyif_t( n_t );

// calculate the maximum density of each species across all cells - move to normal space and try this

// loop through the array to move it to normal space
array<Type> d_gt_csyif(ln_d_gt_csyif.rows(), ln_d_gt_csyif.cols());  // Copy shape
for (int i = 0; i < ln_d_gt_csyif.rows(); i++) {
  for (int j = 0; j < ln_d_gt_csyif.cols(); j++) {
    d_gt_csyif(i,j) = exp(ln_d_gt_csyif(i,j));
  }
}

array<Type> d_gt_rf(ln_d_gt_rf.rows(), ln_d_gt_rf.cols());  // Copy shape
for (int i = 0; i < ln_d_gt_rf.rows(); i++) {
  for (int j = 0; j < ln_d_gt_rf.cols(); j++) {
    d_gt_rf(i,j) = exp(ln_d_gt_rf(i,j));
  }
}

// now get it into a vector
vector<Type> d_gt_csyif_vec = d_gt_csyif.vec();
vector<Type> d_gt_rf_vec = d_gt_rf.vec();

// now get the max value


Type max_csyif = *std::max_element(d_gt_csyif_vec.begin(), d_gt_csyif_vec.end());
Type max_rf = *std::max_element(d_gt_rf_vec.begin(), d_gt_rf_vec.end());

// populate the empty array
for( int t=0; t<n_t; t++){
  for( int g=0; g<n_g_csyif; g++){
    ov_gt_rf_csyif(g,t) = d_gt_rf(g,t) * d_gt_csyif(g,t); // calculate numerator for each cell
    scaled_dens_gt_csyif(g,t) = d_gt_csyif(g,t); // calculate the denominator for each cell
  }
  
  // calculate the index for each year
  ov_rf_csyif_t(t) = ov_gt_rf_csyif.col(t).sum()/scaled_dens_gt_csyif.col(t).sum();
}

// Reporting: overlap
REPORT(ov_rf_csyif_t); // report overlap for full JSOES domain
ADREPORT(ov_rf_csyif_t); // report overlap for full JSOES domain

// REPORT(ov_rf_csyif_t); // report overlap for only coverage overlap between JSOES and PRS
// ADREPORT(ov_rf_csyif_t); // report overlap for only coverage overlap between JSOES and PRS

  
  return jnll;
}