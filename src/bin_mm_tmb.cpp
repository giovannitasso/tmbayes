#include <TMB.hpp>
#include <vector>
// La línea problemática que incluía "density.hpp" se ha eliminado.

template<class Type>
Type objective_function<Type>::operator() () {
  
//==========================
// DATA SECTION
//==========================
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_INTEGER(n_groups); // Number of groups

  
//==========================
// PARAMETER SECTION
//==========================
  PARAMETER_VECTOR(betas);  // fixed effects
  PARAMETER_VECTOR(u);      // random effects
  
  // Parameters for unstructured covariance matrix
  // Assumes 2 random effects per group (intercept and slope)
  PARAMETER_VECTOR(log_stdevs);   // log(stdev) vector [log(sigma_int), log(sigma_slope)]
  PARAMETER(transf_corr);      // Transformed correlation parameter


//==========================
// PRELIMINARY CALCULATIONS
//==========================
  
  // Setup for MVN density
  int n_reff_per_group = 2; 
  density::UNSTRUCTURED_CORR_t<Type> nldens = density::UNSTRUCTURED_CORR_t<Type>(log_stdevs, transf_corr);
  
  // Linear predictor
  vector<Type> eta = X * betas + Z * u;
  vector<Type> p = 1.0 / (1.0 + exp(-eta));

//==========================
// LIKELIHOOD SECTION
//==========================
  Type nll = 0.0;

  // Prior for random effects: u ~ MVN(0, Sigma)
  // We assume u is ordered as [int_g1, ..., int_gN, slope_g1, ..., slope_gN]
  for(int i = 0; i < n_groups; ++i){
    vector<Type> u_group_i(n_reff_per_group);
    u_group_i(0) = u(i);              // Intercept for group i
    u_group_i(1) = u(i + n_groups);   // Slope for group i
    nll += nldens(u_group_i); // Adds *negative* log-density
  }

  // Likelihood for the response: y_i ~ Binomial(1, p_i)
  nll -= sum(dbinom(y, Type(1.0), p, true));

//============================   
//     REPORT section
//============================
  REPORT(betas);
  REPORT(p);       
  REPORT(eta);
  REPORT(u);
  
  // Report variance components
  vector<Type> stdevs = exp(log_stdevs);
  Type rho = nldens.corr(0, 1);
  
  REPORT(stdevs);
  REPORT(rho);
  
  ADREPORT(stdevs);
  ADREPORT(rho);

  return nll;
}
