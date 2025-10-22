#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  
//==========================
// DATA SECTION
//==========================
  DATA_VECTOR(y);          // binary response
  DATA_MATRIX(X);          // design matrix for fixed effects
  DATA_MATRIX(Z);          // random effects design matrix
  // DATA_INTEGER(n_groups); // No necesario para esta versión simple

  
//==========================
// PARAMETER SECTION
//==========================
  PARAMETER_VECTOR(betas);  // fixed effects
  PARAMETER_VECTOR(u);      // random effects
  PARAMETER(logsigma_u);    // log standard deviation of random effects (versión simple)
  // PARAMETER_VECTOR(log_stdevs); // Eliminado temporalmente
  // PARAMETER(transf_corr);      // Eliminado temporalmente


//==========================
// PRELIMINARY CALCULATIONS
//==========================
  // Versión simple:
  Type sigma_u = exp(logsigma_u);
  
  // Linear predictor
  vector<Type> eta = X * betas + Z * u;
  
  // Apply the inverse-logit link function (Cast necesario)
  vector<Type> p = vector<Type>(1.0 / (1.0 + exp(-eta))); 

//==========================
// LIKELIHOOD SECTION
//==========================
  Type nll = 0.0;

  // Prior for random effects (versión simple: dnorm)
  nll -= sum(dnorm(u, Type(0.0), sigma_u, true));

  // Likelihood for the response (sin cambios)
  nll -= sum(dbinom(y, Type(1.0), p, true));

//============================   
//     REPORT section
//============================
  REPORT(betas);
  REPORT(p);       
  REPORT(eta);
  REPORT(u);
  REPORT(sigma_u); // Reportar sigma_u simple
  
  ADREPORT(sigma_u); // ADReport para sigma_u simple

  return nll;
}
