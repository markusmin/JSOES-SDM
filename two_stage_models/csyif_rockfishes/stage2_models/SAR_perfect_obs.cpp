
#include <TMB.hpp>
#include <algorithm>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{

using namespace density;
  // SAR:
  // Data
  DATA_VECTOR(z_i);  // survival of fish i
  DATA_IVECTOR(t_i); // index for the year of fish i
  DATA_INTEGER(n_t); // number of years in the dataset
  
  DATA_VECTOR(ov_ab_t);  // overlap in year t
  // DATA_MATRIX(Sigma_ov_ab);  // covariance matrix for overlap
  
  // Parameters
  // PARAMETER(ln_sigma_t); // standard deviation for random effect of year
  PARAMETER(beta_0); // this is the intercept
  PARAMETER(beta_ov_ab); // this is the effect of biomass-weighted overlap
  
  
  // Random effects
  // PARAMETER_VECTOR(epsilon_t); // random effect of year
  
  // PARAMETER_VECTOR(ov_ab_t_latent); // actual overlap
  
  // Objective function
  Type jnll = 0;
  
  // Estimate actual overlap as a latent random effect
  // MVN for errors in variables of overlap metric
  // vector<Type> residual(n_t);
  // MVNORM_t<Type> neg_log_dmvnorm(Sigma_ov_ab);
  // residual = vector<Type>(ov_ab_t - ov_ab_t_latent); // here we are calculating the residuals between the data and the underlying latent state
  // jnll += neg_log_dmvnorm(residual); // then we use those residuals and evaluate them
  
  // for( int t=0; t<n_t; t++){
  //   // first calculate the residuals - as this is how we can use MVNORM_t which assumes mean of zero
  //   residual = vector<Type>(ov_ab_t(t));
  //   
  //   jnll += neg_log_dmvnorm();
  // }
  
  
  // for( int t=0; t<n_t; t++){
  //   // jnll -= dnorm(ov_ab_t(t), ov_ab_t_latent(t), sigma_ov_ab_t(t));
  //   jnll += MVNORM_t(Sigma_ov_ab)(ov_ab_t(t));
  // }
  // MVNORM_t<Type> neg_log_density(Sigma_ov_ab);
  // for( int t=0; t<n_t; t++){
  //   jnll += neg_log_density(ov_ab_t(t));
  // }
  // jnll += MVNORM(Sigma_ov_ab)(ov_ab_t);
  
  
  // SAR: random effect of year
  // jnll -= sum(dnorm(epsilon_t, 0, exp(ln_sigma_t)));
  
  
  // SAR model, using biomass-weighted overlap index
  for( int i=0; i<z_i.size(); i++){
    // jnll -= dbinom_robust(z_i(i), Type(1), beta_ov_ab*ov_ab_t_latent(t_i(i)) + epsilon_t(t_i(i)), true ); // this is a Bernoulli (by using Type(1) as the size argument)
    jnll -= dbinom_robust(z_i(i), Type(1), beta_0 + beta_ov_ab*ov_ab_t(t_i(i)), true ); // this is a Bernoulli (by using Type(1) as the size argument)
  }
  
  
  // Reporting: SAR
  // REPORT( ln_sigma_t );
  // REPORT( epsilon_t );
  REPORT( beta_ov_ab );
  REPORT( beta_0 );

  
  return jnll;
}