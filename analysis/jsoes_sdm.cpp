#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  // Data
  DATA_VECTOR( N_i );  // density for observation i
  DATA_VECTOR( W_i );  // weight for observation i
  DATA_VECTOR( simulate_t ); // vector of zeros or 1s for which to simulate new values
  
  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  
  // Projection matrices
  DATA_SPARSE_MATRIX(A_is);
  
  // Parameters
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(ln_tau);
  PARAMETER(ln_kappa);
  // PARAMETER(ln_sigma);
  
  // Random effects
  PARAMETER_VECTOR(omega_s);
  PARAMETER_VECTOR(psi_s);
  
  // Objective function
  Type jnll = 0;
  
  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));
  
  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = (exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2) * exp(2*ln_tau);
  jnll += SCALE( GMRF(Q), 1/exp(ln_tau) )( omega_s );
  jnll += SCALE( GMRF(Q), 1/exp(ln_tau) )( psi_s );
  
  SIMULATE{
    // simulate omega
    SCALE( GMRF(Q), 1/exp(ln_tau) ).simulate(omega_s);
    // omega_s = GMRF(Q).simulate();
    // simulate psi
    SCALE( GMRF(Q), 1/exp(ln_tau) ).simulate(psi_s);
    // psi_s = GMRF(Q).simulate();
  }
  // Project using bilinear interpolation
  vector<Type> omega_i( A_is.rows() );
  omega_i = A_is * omega_s;
  vector<Type> psi_i( A_is.rows() );
  psi_i = A_is * psi_s;
  
  // store lambda_s for later calculations
  vector<Type> log_lambda_s( A_is.rows() );
  for( int i=0; i<N_i.size(); i++){
    log_lambda_s(i) = alpha + omega_s(i);
  }
  // create intermediate variable log_lambda_i (map log_lambda_s onto samples
  vector<Type> log_lambda_i( A_is.rows() );
  
  // Probability of data conditional on random effects
  // counts
  for( int i=0; i<N_i.size(); i++){
    if( !R_IsNA(asDouble(N_i(i))) ){
      log_lambda_i(i) = alpha + omega_i(i);
      jnll -= dpois( N_i(i), exp(log_lambda_i(i)), true );
    }
  }
  
  // weights
  // create a vector to store average body size in each cell
  vector<Type> log_mu_s( A_is.rows() );
  for( int i=0; i<W_i.size(); i++){
    log_mu_s(i) = beta + 0.5*omega_s(i) + 0.5*psi_s(i);
  }
  
  // create intermediate variable log_mu_i (map log_mu_s onto samples)
  vector<Type> log_mu_i( A_is.rows() );
  
  // set sigma to 0.5
  Type sigma = 0.5;
  
  // evaluate probability of weights
  for( int i=0; i<W_i.size(); i++){
    if( !R_IsNA(asDouble(W_i(i))) ){
      log_mu_i(i) = beta + 0.5*omega_i(i) + 0.5*psi_i(i);
      jnll -= dgamma( W_i(i), pow(sigma, -2), exp(log_mu_i(i))*pow(sigma, 2), true );
    }
  }
  
  
  // simulate sampling design
  for( int t=0; t<N_i.size(); t++){
    if( simulate_t(t) == 1 ){
      SIMULATE{
        N_i(t) = rpois(exp(log_lambda_s(t)));
        W_i(t) = rgamma(pow(sigma, -2), exp(log_mu_s(t))*pow(sigma, 2));
      }
    }
  }
  
  
  // W_bar 3: Area-weighted model for body size
  Type W_bar3 = 0;
  // Plug-in estimator
  W_bar3 = exp(log_mu_s).sum()/100;
  
  
  // W_bar 4: Abundance-weighted model for body size
  Type W_bar4 = 0;
  // Plug-in estimator
  W_bar4 = (exp(log_lambda_s)/exp(log_lambda_s).sum()*exp(log_mu_s)).sum();
  
  
  // Report
  REPORT( Q );
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( log_lambda_s );
  REPORT( log_mu_s );
  REPORT( omega_s );
  REPORT( psi_s );
  // Use ADREPORT to get automated bias correction for model-based estimators
  REPORT( W_bar3 );
  REPORT( W_bar4 );
  ADREPORT( W_bar3 );
  ADREPORT( W_bar4 );
  
  
  // Report simulated quantities
  SIMULATE{ REPORT( log_lambda_s ); } // report numerical density in each cell (log_lambda_s)
  SIMULATE{ REPORT( log_mu_s ); } // report average body size in each cell (log_mu_s)
  SIMULATE{ REPORT( N_i ); } // report simulated counts
  SIMULATE{ REPORT( W_i ); } // report simulated weights
  SIMULATE{ REPORT( omega_s ); } // report simulated omega_s
  SIMULATE{ REPORT( psi_s ); } // report simulated psi_s
  
  // Reporting
  return jnll;
  
}




