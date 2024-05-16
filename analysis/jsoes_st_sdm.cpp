
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
  
  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  
  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);
  
  // Parameters
  PARAMETER_VECTOR( beta_t );
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );
  PARAMETER( ln_phi ); // phi term in tweedie
  PARAMETER( finv_power ); // power parameter in tweedie
  
  // Random effects
  PARAMETER_VECTOR( omega_s );
  PARAMETER_VECTOR( epsilon_st );
  
  // Objective funcction
  Type jnll = 0;
  
  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));
  
  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  jnll += SCALE( GMRF(Q), 1/exp(ln_tau) )( omega_s );
  
  // Probability of data conditional on random effects
  vector<Type> omega_i( A_is.rows() );
  omega_i = A_is * omega_s;
  
  vector<Type> omega_g( A_gs.rows() );
  omega_g = A_gs * omega_s;
  
  int n_i = A_is.rows();
  int n_g = A_gs.rows();
  vector<Type> dhat_i( n_i );
  dhat_i.setZero();
  
  
  for( int i=0; i<D_i.size(); i++){
    dhat_i(i) = exp( beta_t(t_i(i)) + omega_i(i));
    jnll -= dtweedie( D_i(i), dhat_i(i), exp(ln_phi), Type(1.0)+invlogit(finv_power), true );
  }
  
  // modeled density in the projection grid
  array<Type> ln_d_gt( n_g, n_t );
  for( int t=0; t<n_t; t++){
    for( int g=0; g<n_g; g++){
      ln_d_gt(g,t) = beta_t(t) + omega_g(g);
    }
  }
  
  
  // Reporting
  REPORT( dhat_i );
  REPORT( beta_t );
  REPORT( ln_tau );
  REPORT( ln_kappa );
  REPORT( ln_phi );
  REPORT( finv_power );
  REPORT(SigmaE);
  REPORT(Range);
  REPORT(ln_d_gt);
  
  return jnll;
}