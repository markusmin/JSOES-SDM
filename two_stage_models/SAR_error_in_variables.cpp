
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
  
  DATA_VECTOR(o_t);  // overlap in year t
  DATA_VECTOR(o_t);  // uncertainty in oerlap
  
  // Parameters
  PARAMETER(ln_sigma_t);
  PARAMETER(beta_ov_ab); // this is the effect of biomass-weighted overlap
  
  // Random effects
  PARAMETER_VECTOR(epsilon_t); // random effect of year
  
  
  // SAR model, using biomass-weighted overlap index
  for( int i=0; i<z_i.size(); i++){
    jnll -= dbinom_robust(z_i(i), Type(1), beta_ov_ab*ov_ab_t(t_i(i)) + epsilon_t(t_i(i)), true ); // this is a Bernoulli (by using Type(1) as the size argument)
  }
  
  
  // Reporting: SAR
  REPORT( ln_sigma_t );
  REPORT( epsilon_t );
  REPORT( beta_ov_ab );
  
  
  

  
  return jnll;
}