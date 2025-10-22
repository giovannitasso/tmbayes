// Keep includes minimal: TMB.hpp should handle dependencies.
#include <TMB.hpp>

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
  PARAMETER(transf_corr);      // Transformed correlation parameter (-inf, +inf)


//==========================
// PRELIMINARY CALCULATIONS
//==========================
  
  // Setup for MVN density using TMB's built-in tools
  int n_reff_per_group = 2; // Hardcoded for intercept and slope for now
  density::UNSTRUCTURED_CORR_t<Type> nldens = density::UNSTRUCTURED_CORR_t<Type>(log_stdevs, transf_corr);
  
  // Linear predictor (using TMB/Eigen matrix operations)
  vector<Type> eta = X * betas + Z * u;
  
  // Apply the inverse-logit link function to get probabilities
  // Explicit cast to vector<Type> to avoid Eigen lazy-evaluation issues
  vector<Type> p = vector<Type>(1.0 / (1.0 + exp(-eta))); 

//==========================
// LIKELIHOOD SECTION
//==========================
  Type nll = 0.0; // Initialize negative log-likelihood

  // Prior for random effects: u ~ MVN(0, Sigma)
  // Loop through each group and apply the MVN density
  // Assumes u is ordered as [int_g1, ..., int_gN, slope_g1, ..., slope_gN]
  for(int i = 0; i < n_groups; ++i){
    vector<Type> u_group_i(n_reff_per_group);
    u_group_i(0) = u(i);              // Intercept for group i
    u_group_i(1) = u(i + n_groups);   // Slope for group i
    nll += nldens(u_group_i); // Adds the *negative* log-density from the MVN prior
  }

  // Likelihood for the response: y_i ~ Binomial(1, p_i)
  // Use TMB's vectorized dbinom, subtracting log-prob because it's NLL
  nll -= sum(dbinom(y, Type(1.0), p, true));

//============================   
//     REPORT section
//============================
  REPORT(betas);
  REPORT(p);       
  REPORT(eta);
  REPORT(u);
  
  // Report the standard deviations and correlation derived from parameters
  vector<Type> stdevs = exp(log_stdevs);
  Type rho = nldens.corr(0, 1); // Extract correlation between effect 1 and 2
  
  REPORT(stdevs);
  REPORT(rho);
  
  // Calculate standard errors for derived quantities using the delta method
  ADREPORT(stdevs);
  ADREPORT(rho);

  return nll;
}
