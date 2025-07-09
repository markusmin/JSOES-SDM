
#include <TMB.hpp>
#include <algorithm>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{

using namespace density;

  // Data
  DATA_VECTOR(D_j_a);  // density of species A for measurement j
  DATA_VECTOR(D_k_b);  // density of species B for measurement k
  DATA_IVECTOR(t_j); // index for the year of measurement j
  DATA_IVECTOR(t_k); // index for the year of measurement k
  DATA_INTEGER(n_t); // number of years in the dataset
  DATA_VECTOR(temp_j);  // temperature for measurement j
  DATA_VECTOR(temp_k);  // temperature for measurement k
  DATA_VECTOR(dist_j);  // distance from shore for measurement j
  DATA_VECTOR(dist_k);  // distance from shore for measurement k
  
  // SPDE objects for species A
  DATA_SPARSE_MATRIX(M0_a);
  DATA_SPARSE_MATRIX(M1_a);
  DATA_SPARSE_MATRIX(M2_a);
  
  // SPDE objects for species B
  DATA_SPARSE_MATRIX(M0_b);
  DATA_SPARSE_MATRIX(M1_b);
  DATA_SPARSE_MATRIX(M2_b);
  
  // Projection matrices for species A
  // A_is is the projection matrix from vertices to samples, and therefore has unique dimensions for each survey
  DATA_SPARSE_MATRIX(A_is);
  // A_gs is the projection matrix from the vertices to the projection grid for the survey domain, where the grid dimensions are the same across surveys/species
  DATA_SPARSE_MATRIX(A_gs);
  
  // Projection matrices for species B
  // B_is is the projection matrix from vertices to samples, and therefore has unique dimensions for each survey
  DATA_SPARSE_MATRIX(B_is);
  // B_gs is the projection matrix from the vertices to the projection grid for the survey domain, where the grid dimensions are the same across surveys/species
  DATA_SPARSE_MATRIX(B_gs);
  
  DATA_MATRIX(temp_gt); // temperature at each location in each year
  DATA_VECTOR(dist_g); // distance from shore at each location
  
  
  // Parameters
  
  // Parameters for SDM for species A
  PARAMETER_VECTOR( beta_t_a );
  PARAMETER( beta_temp_a );
  PARAMETER( beta_dist_a );
  PARAMETER( ln_tau_omega_a );
  PARAMETER( ln_tau_epsilon_a );
  PARAMETER( ln_kappa_a );
  PARAMETER( ln_phi_a ); // phi term in tweedie
  PARAMETER( finv_power_a ); // power parameter in tweedie
  
  // Parameters for SDM for species B
  PARAMETER_VECTOR( beta_t_b );
  PARAMETER( beta_temp_b );
  PARAMETER( beta_dist_b );
  PARAMETER( ln_tau_omega_b );
  PARAMETER( ln_tau_epsilon_b );
  PARAMETER( ln_kappa_b );
  PARAMETER( ln_phi_b ); // phi term in tweedie
  PARAMETER( finv_power_b ); // power parameter in tweedie
  
  // Random effects
  
  // for species A
  PARAMETER_VECTOR( omega_s_a );
  PARAMETER_MATRIX( epsilon_st_a );
  
  // for species B
  PARAMETER_VECTOR( omega_s_b );
  PARAMETER_MATRIX( epsilon_st_b );
  
  
  // Objective function
  Type jnll = 0;
  int n_i_a = A_is.rows();
  int n_g_a = A_gs.rows();
  int n_i_b = B_is.rows();
  int n_g_b = B_gs.rows();
  
  // Derived quantities
  // SDM derived quantities for species A
  Type Range_a = sqrt(8) / exp( ln_kappa_a );
  Type SigmaO_a = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_a) * exp(2*ln_kappa_a));
  Type SigmaE_a = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_a) * exp(2*ln_kappa_a));
  
  // SDM derived quantities for species B
  Type Range_b = sqrt(8) / exp( ln_kappa_b );
  Type SigmaO_b = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega_b) * exp(2*ln_kappa_b));
  Type SigmaE_b = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon_b) * exp(2*ln_kappa_b));
  
  // Probability of random effects
  
  // SDM: species A
  Eigen::SparseMatrix<Type> Q_a = exp(4*ln_kappa_a)*M0_a + Type(2.0)*exp(2*ln_kappa_a)*M1_a + M2_a;
  // spatial random effect - scaled by ln_tau_omega
  jnll += SCALE( GMRF(Q_a), 1/exp(ln_tau_omega_a) )( omega_s_a );
  // spatio-temporal random effect - scaled by ln_tau_epsilon
  for( int t=0; t<n_t; t++){
    jnll += SCALE( GMRF(Q_a), 1/exp(ln_tau_epsilon_a) )( epsilon_st_a.col(t) );
  }
  
  // SDM: species B
  Eigen::SparseMatrix<Type> Q_b = exp(4*ln_kappa_b)*M0_b + Type(2.0)*exp(2*ln_kappa_b)*M1_b + M2_b;
  // spatial random effect - scaled by ln_tau_omega
  jnll += SCALE( GMRF(Q_b), 1/exp(ln_tau_omega_b) )( omega_s_b );
  // spatio-temporal random effect - scaled by ln_tau_epsilon
  for( int t=0; t<n_t; t++){
    jnll += SCALE( GMRF(Q_b), 1/exp(ln_tau_epsilon_b) )( epsilon_st_b.col(t) );
  }
  
  // SDM projections for species A
  vector<Type> omega_g_a( n_g_a );
  omega_g_a = A_gs * omega_s_a;
  matrix<Type> epsilon_gt_a( n_g_a, n_t );
  epsilon_gt_a = A_gs * epsilon_st_a;
  
  // SDM projections for species B
  vector<Type> omega_g_b( n_g_b );
  omega_g_b = B_gs * omega_s_b;
  matrix<Type> epsilon_gt_b( n_g_b, n_t );
  epsilon_gt_b = B_gs * epsilon_st_b;
  
  // Probability of data conditional on random effects
  
  // SDM species A
  vector<Type> omega_i_a( n_i_a );
  omega_i_a = A_is * omega_s_a;
  matrix<Type> epsilon_it_a( n_i_a, n_t );
  epsilon_it_a = A_is * epsilon_st_a;

  vector<Type> dhat_j_a( n_i_a );
  dhat_j_a.setZero();


  for( int j=0; j<D_j_a.size(); j++){
    dhat_j_a(j) = exp( beta_t_a(t_j(j)) + beta_temp_a*temp_j(j) + beta_dist_a*dist_j(j) + omega_i_a(j) + epsilon_it_a(j,t_j(j)) );
    jnll -= dtweedie( D_j_a(j), dhat_j_a(j), exp(ln_phi_a), Type(1.0)+invlogit(finv_power_a), true );
  }
  
  // SDM species B
  vector<Type> omega_i_b( n_i_b );
  omega_i_b = B_is * omega_s_b;
  matrix<Type> epsilon_it_b( n_i_b, n_t );
  epsilon_it_b = B_is * epsilon_st_b;

  vector<Type> dhat_k_b( n_i_b );
  dhat_k_b.setZero();


  for( int k=0; k<D_k_b.size(); k++){
    dhat_k_b(k) = exp( beta_t_b(t_k(k)) + beta_temp_b*temp_k(k) + beta_dist_b*dist_k(k) + omega_i_b(k) + epsilon_it_b(k,t_k(k)) );
    jnll -= dtweedie( D_k_b(k), dhat_k_b(k), exp(ln_phi_b), Type(1.0)+invlogit(finv_power_b), true );
  }
  
  
  
  // SDM: modeled density in the projection grid for species A
  // Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
  array<Type> ln_d_gt_a( n_g_a, n_t );
  for( int t=0; t<n_t; t++){
    for( int g=0; g<n_g_a; g++){
      ln_d_gt_a(g,t) = beta_t_a(t) + beta_temp_a*temp_gt(g,t) + beta_dist_a*dist_g(g) + omega_g_a(g) + epsilon_gt_a(g,t);
    }
  }
  
  // SDM: modeled density in the projection grid for species B
  // Note that the projection grid is the same for all species - so temp_gt and dist_g are universal, and the number of years (t) are the same
  array<Type> ln_d_gt_b( n_g_b, n_t );
  for( int t=0; t<n_t; t++){
    for( int g=0; g<n_g_b; g++){
      ln_d_gt_b(g,t) = beta_t_b(t) + beta_temp_b*temp_gt(g,t) + beta_dist_b*dist_g(g) + omega_g_b(g) + epsilon_gt_b(g,t);
    }
  }
  
  
  // Use the modeled density in the projection grid to calculate biomass-weighted overlap
  // Use a nested for loop to loop through all years and all grid cells
  // note that n_g_a and n_g_b should have the same dimensions and therefore we can use either to construct the array
  
  // construct two arrays: one to store the numerator and one to store the denominator for the index
  array<Type> ov_gt_ab( n_g_a, n_t );
  array<Type> scaled_dens_gt_b( n_g_a, n_t );
  
  // construct a matrix to store the overlap index in each year
  vector<Type> ov_ab_t( n_t );
  
  // calculate the maximum density of each species across all cells - move to normal space and try this
  
  // loop through the array to move it to normal space
  array<Type> d_gt_a(ln_d_gt_a.rows(), ln_d_gt_a.cols());  // Copy shape
  for (int i = 0; i < ln_d_gt_a.rows(); i++) {
    for (int j = 0; j < ln_d_gt_a.cols(); j++) {
      d_gt_a(i,j) = exp(ln_d_gt_a(i,j));  
    }
  }
  
  array<Type> d_gt_b(ln_d_gt_b.rows(), ln_d_gt_b.cols());  // Copy shape
  for (int i = 0; i < ln_d_gt_b.rows(); i++) {
    for (int j = 0; j < ln_d_gt_b.cols(); j++) {
      d_gt_b(i,j) = exp(ln_d_gt_b(i,j));  
    }
  }
  
  // now get it into a vector
  vector<Type> d_gt_a_vec = d_gt_a.vec();
  vector<Type> d_gt_b_vec = d_gt_b.vec();
  
  // now get the max value
  
  Type max_A = *std::max_element(d_gt_a_vec.begin(), d_gt_a_vec.end());
  Type max_B = *std::max_element(d_gt_b_vec.begin(), d_gt_b_vec.end());

  

  // populate the empty array
  for( int t=0; t<n_t; t++){
    for( int g=0; g<n_g_a; g++){
      ov_gt_ab(g,t) = d_gt_a(g,t)/max_A * d_gt_b(g,t)/max_B; // calculate numerator for each cell
      scaled_dens_gt_b(g,t) = d_gt_b(g,t)/max_B; // calculate the denominator for each cell
    }
    
    // calculate the index for each year
    ov_ab_t(t) = ov_gt_ab.col(t).sum()/scaled_dens_gt_b.col(t).sum();
  }
  
  
  // Reporting: SDM
  REPORT( dhat_j_a );
  REPORT( dhat_k_b );
  REPORT( beta_t_a );
  REPORT( beta_t_b );
  REPORT( beta_temp_a );
  REPORT( beta_temp_b );
  REPORT( beta_dist_a );
  REPORT( beta_dist_b );
  REPORT( ln_tau_omega_a );
  REPORT( ln_tau_omega_b );
  REPORT( ln_tau_epsilon_a );
  REPORT( ln_tau_epsilon_b );
  REPORT( ln_kappa_a );
  REPORT( ln_kappa_b );
  REPORT( ln_phi_a );
  REPORT( ln_phi_b );
  REPORT( finv_power_a );
  REPORT( finv_power_b );
  REPORT(SigmaO_a);
  REPORT(SigmaE_a);
  REPORT(Range_a);
  REPORT(SigmaO_b);
  REPORT(SigmaE_b);
  REPORT(Range_b);
  REPORT(ln_d_gt_a);
  REPORT(ln_d_gt_b);
  
  // report the overlap by year
  REPORT(ov_ab_t);
  
  return jnll;
}