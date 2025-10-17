#' @title Fit a Binomial Generalized Linear Mixed Model with TMB
#' @description This function fits a binomial GLMM using a pre-compiled TMB model.
#' It can handle both random intercepts and random slopes, depending on the
#' structure of the design matrix Z.
#'
#' @param y A numeric vector representing the binary response variable (0s and 1s).
#' @param X A numeric matrix of fixed effects covariates. An intercept is recommended.
#' @param Z A numeric matrix for the random effects design. For a simple random
#'        intercept model, this is the model matrix for the group intercepts. For a
#'        random slope model, this matrix will have more columns.
#' @param initial_betas A numeric vector for the starting values of the fixed effects (betas).
#'        If NULL, defaults to zeros.
#' @param initial_logsigma_u A single numeric value for the starting value of the log standard
#'        deviation of the random effects. Defaults to 0.
#'
#' @return A list object of class `tmbayes_fit` containing the model results.
#' @export
#' @examples
#' \dontrun{
#' # --- Example 1: Random Intercepts Only ---
#' library(tmbayes)
#' set.seed(1234)
#' n_groups <- 10
#' n_per_group <- 20
#' n_tot <- n_groups * n_per_group
#' 
#' group <- rep(1:n_groups, each = n_per_group)
#' x <- rnorm(n_tot)
#' u_true <- rnorm(n_groups, 0, 1.0) # Random intercepts
#' 
#' eta <- -1.0 + 1.5 * x + u_true[group]
#' p <- 1 / (1 + exp(-eta))
#' y <- rbinom(n_tot, size = 1, prob = p)
#' 
#' X <- cbind(1, x)
#' Z_intercept <- model.matrix(~ factor(group) - 1)
#' 
#' fit_intercepts <- fit_binomial_glmm(y = y, X = X, Z = Z_intercept)
#' print(fit_intercepts)
#'
#' # --- Example 2: Random Intercepts and Random Slopes ---
#' set.seed(1234)
#' # Simulate random slopes for x that vary by group
#' u_slope_true <- rnorm(n_groups, 0, 0.5)
#' 
#' # Linear predictor now includes a random slope component
#' eta_slopes <- -1.0 + 1.5 * x + u_true[group] + u_slope_true[group] * x
#' p_slopes <- 1 / (1 + exp(-eta_slopes))
#' y_slopes <- rbinom(n_tot, size = 1, prob = p_slopes)
#'
#' # Build the Z matrix for random intercepts AND slopes
#' # We need to construct this matrix manually.
#' # Each observation needs its group's intercept and its group's slope effect.
#' Z_slopes <- matrix(0, nrow = n_tot, ncol = n_groups * 2)
#' for (i in 1:n_tot) {
#'   group_idx <- group[i]
#'   # Column for the intercept of the i-th observation's group
#'   Z_slopes[i, group_idx] <- 1
#'   # Column for the slope of the i-th observation's group
#'   Z_slopes[i, n_groups + group_idx] <- x[i]
#' }
#' 
#' # Fit the model
#' # The C++ code assumes all random effects (intercepts and slopes)
#' # share the same standard deviation sigma_u.
#' fit_slopes <- fit_binomial_glmm(y = y_slopes, X = X, Z = Z_slopes)
#' print(fit_slopes)
#' }
fit_binomial_glmm <- function(y, X, Z, initial_betas = NULL, initial_logsigma_u = 0) {
  # ... (el resto de tu función no necesita cambios)
  # ---- 1. Data Validation and Preparation ----
  if (!is.numeric(y) || !all(y %in% c(0, 1))) {
    stop("Response variable 'y' must be a numeric vector of 0s and 1s.")
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("'X' must be a numeric matrix for fixed effects.")
  }
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("'Z' must be a numeric matrix for random effects.")
  }
  if (length(y) != nrow(X) || length(y) != nrow(Z)) {
    stop("Number of rows in 'y', 'X', and 'Z' must be equal.")
  }
  
  tmb_data <- list(
    y = y, 
    X = X, 
    Z = Z
  )
  
  # ---- 2. Parameter Initialization ----
  if (is.null(initial_betas)) {
    initial_betas <- rep(0, ncol(X))
  }
  
  tmb_params <- list(
    betas = initial_betas, 
    u = rep(0, ncol(Z)),
    logsigma_u = initial_logsigma_u
  )
  
  # ---- 3. TMB Model Fitting ----
  obj <- TMB::MakeADFun(
    data = tmb_data, 
    parameters = tmb_params,
    random = "u",
    DLL = "tmbayes",
    silent = TRUE
  )
  
  opt <- stats::nlminb(
    start = obj$par, 
    objective = obj$fn, 
    gradient = obj$gr
  )
  
  # ---- 4. Result Extraction and Formatting ----
  report <- TMB::sdreport(obj)
  summary_report <- summary(report)
  
  num_fixed_effects <- ncol(X)
  fixed_eff_idx <- which(rownames(summary_report) == "betas")
  
  correct_fixed_eff_idx <- fixed_eff_idx[1:num_fixed_effects]
  fixed_coeffs <- summary_report[correct_fixed_eff_idx, , drop = FALSE]
  
  rownames(fixed_coeffs) <- if (!is.null(colnames(X))) colnames(X) else paste0("beta_", 1:nrow(fixed_coeffs))
  
  random_eff_idx <- which(rownames(summary_report) == "sigma_u")
  random_coeffs <- summary_report[random_eff_idx, , drop = FALSE]
  rownames(random_coeffs) <- "sigma_u"
  
  final_summary <- rbind(fixed_coeffs, random_coeffs)
  colnames(final_summary) <- c("Estimate", "Std. Error")
  
  final_summary <- as.data.frame(final_summary)
  
  final_summary$"z value" <- final_summary[, "Estimate"] / final_summary[, "Std. Error"]
  final_summary$"Pr(>|z|)" <- 2 * stats::pnorm(-abs(final_summary$"z value"))
  
  output <- list(
    summary = final_summary,
    report = report,
    opt = opt,
    obj = obj
  )
  
  class(output) <- "tmbayes_fit" 
  
  return(output)
}


#' @title Print method for tmbayes_fit objects
#' @description A custom print method to display a concise summary of the fitted model.
#' @param x An object of class `tmbayes_fit`.
#' @param ... Additional arguments passed to `print`.
#' @method print tmbayes_fit
#' @export
print.tmbayes_fit <- function(x, ...) {

  cat("Laplace Approximation Fit from 'tmbayes'\n\n") 
  

  cat("Formula: Binomial GLMM\n")
  cat("Optimizer status:", x$opt$message, "\n\n")
  
  cat("Parameter Estimates:\n")
  print(round(x$summary, 4))
  cat("\n")
}
