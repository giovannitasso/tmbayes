#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  
//==========================
// DATA SECTION
//==========================
  DATA_VECTOR(y);          // binary response
  DATA_MATRIX(X);          // design matrix for fixed effects
  DATA_MATRIX(Z);          // random effects design matrix

  
//==========================
// PARAMETER SECTION
//==========================
  PARAMETER_VECTOR(betas);  // fixed effects
  PARAMETER_VECTOR(u);      // random effects
  PARAMETER(logsigma_u);    // log standard deviation of random effects

//==========================
// PRELIMINARY CALCULATIONS
//==========================
  // Transform log-standard deviation to standard deviation
  Type sigma_u = exp(logsigma_u);
  
  // Calculate the linear predictor
  // eta = X*betas + Z*u
  vector<Type> eta = X * betas + Z * u;
  
  // Apply the inverse-logit link function to get probabilities
  vector<Type> p = 1.0 / (1.0 + exp(-eta));

//==========================
// LIKELIHOOD SECTION
//==========================
  // Initialize the negative log-likelihood
  Type nll = 0.0;

  // Prior for random effects: u ~ N(0, sigma_u^2)
  // dnorm(x, mean, sd, give_log)
  // The 'true' at the end means it returns the log-density.
  // We subtract because we are calculating the *negative* log-likelihood.
  nll -= sum(dnorm(u, Type(0.0), sigma_u, true));

  // Likelihood for the response: y_i ~ Binomial(1, p_i)
  // dbinom(x, size, prob, give_log)
  // We use Type(1.0) for the size of the binomial trial.
  nll -= sum(dbinom(y, Type(1.0), p, true));

//============================   
//     REPORT section
//============================
  // These quantities can be extracted from the model object in R
  REPORT(betas);
  REPORT(sigma_u);
  REPORT(p);       
  REPORT(eta);
  REPORT(u);
  
  // ADREPORT will calculate standard errors for these using the delta method
  ADREPORT(sigma_u);

  return nll;
}
