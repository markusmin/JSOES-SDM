
#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data
  DATA_VECTOR(D_i);  // density for measurement i
  DATA_IVECTOR(t_i); // index for the year of measurement i
  DATA_INTEGER(n_t); // number of years in the dataset
  DATA_VECTOR(temp_i);  // temperature for measurement i
  DATA_VECTOR(dist_i);  // distance from shore for measurement i
  
  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  
  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);
  
  DATA_MATRIX(temp_gt); // temperature at each location in each year
  DATA_VECTOR(dist_g); // distance from shore at each location
  
  
  // Parameters
  PARAMETER_VECTOR( beta_t );
  PARAMETER( beta_temp );
  PARAMETER( beta_dist );
  PARAMETER( ln_tau_omega );
  PARAMETER( ln_tau_epsilon );
  PARAMETER( ln_kappa );
  PARAMETER( ln_phi ); // phi term in tweedie
  PARAMETER( finv_power ); // power parameter in tweedie
  
  // Random effects
  PARAMETER_VECTOR( omega_s );
  PARAMETER_MATRIX( epsilon_st );
  
  // Objective funcction
  Type jnll = 0;
  int n_i = A_is.rows();
  int n_g = A_gs.rows();
  
  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_omega) * exp(2*ln_kappa));
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau_epsilon) * exp(2*ln_kappa));
  
  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  // spatial random effect - scaled by ln_tau_omega
  jnll += SCALE( GMRF(Q), 1/exp(ln_tau_omega) )( omega_s );
  // spatio-temporal random effect - scaled by ln_tau_epsilon
  for( int t=0; t<n_t; t++){
    jnll += SCALE( GMRF(Q), 1/exp(ln_tau_epsilon) )( epsilon_st.col(t) );
  }
  
  // projections
  vector<Type> omega_g( n_g );
  omega_g = A_gs * omega_s;
  matrix<Type> epsilon_gt( n_g, n_t );
  epsilon_gt = A_gs * epsilon_st;
  
  // Probability of data conditional on random effects
  vector<Type> omega_i( n_i );
  omega_i = A_is * omega_s;
  matrix<Type> epsilon_it( n_i, n_t );
  epsilon_it = A_is * epsilon_st;
  
  vector<Type> dhat_i( n_i );
  dhat_i.setZero();
  
  
  for( int i=0; i<D_i.size(); i++){
    dhat_i(i) = exp( beta_t(t_i(i)) + beta_temp*temp_i(i) + beta_dist*dist_i(i) + omega_i(i) + epsilon_it(i,t_i(i)) );
    jnll -= dtweedie( D_i(i), dhat_i(i), exp(ln_phi), Type(1.0)+invlogit(finv_power), true );
  }
  
  // modeled density in the projection grid
  array<Type> ln_d_gt( n_g, n_t );
  for( int t=0; t<n_t; t++){
    for( int g=0; g<n_g; g++){
      ln_d_gt(g,t) = beta_t(t) + beta_temp*temp_gt(g,t) + beta_dist*dist_g(g) + omega_g(g) + epsilon_gt(g,t);
    }
  }
  
  
  // Reporting
  REPORT( dhat_i );
  REPORT( beta_t );
  REPORT( beta_temp );
  REPORT( beta_dist );
  REPORT( ln_tau_omega );
  REPORT( ln_tau_epsilon );
  REPORT( ln_kappa );
  REPORT( ln_phi );
  REPORT( finv_power );
  REPORT(SigmaO);
  REPORT(SigmaE);
  REPORT(Range);
  REPORT(ln_d_gt);
  
  return jnll;
}