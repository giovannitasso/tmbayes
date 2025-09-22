#' @title Fit a Binomial Generalized Linear Mixed Model with TMB
#' @description This function provides a user-friendly interface to fit a binomial GLMM
#' using a pre-compiled TMB model. It handles data preparation, parameter initialization,
#' model fitting, and result extraction.
#'
#' @param y A numeric vector representing the binary response variable (0s and 1s).
#' @param X A numeric matrix of fixed effects covariates. An intercept is recommended.
#' @param Z A numeric matrix for the random effects design. Typically a model matrix for group intercepts.
#' @param initial_betas A numeric vector for the starting values of the fixed effects (betas).
#'        If NULL, defaults to zeros.
#' @param initial_logsigma_u A single numeric value for the starting value of the log standard
#'        deviation of the random effects. Defaults to 0.
#'
#' @return A list object of class `tmbayes_fit` containing:
#' \itemize{
#'   \item `summary`: A data frame with estimates, standard errors, and z-values for fixed and random effects parameters.
#'   \item `report`: The full report from `TMB::sdreport`.
#'   \item `opt`: The raw optimization output from `nlminb`.
#'   \item `obj`: The TMB object created by `TMB::MakeADFun`.
#' }
#' @export
#' @examples
#' \dontrun{
#' # 1. Simulate data
#' library(tmbayes) # Note: library name updated
#' set.seed(1234)
#' n_groups <- 10
#' n_per_group <- 20
#' n_tot <- n_groups * n_per_group
#' 
#' group <- rep(1:n_groups, each = n_per_group)
#' x <- rnorm(n_tot)
#' u_true <- rnorm(n_groups, 0, 1.0)
#' 
#' eta <- -1.0 + 1.5 * x + u_true[group]
#' p <- 1 / (1 + exp(-eta))
#' y <- rbinom(n_tot, size = 1, prob = p)
#' 
#' X <- cbind(1, x)
#' Z <- model.matrix(~ factor(group) - 1)
#' 
#' # 2. Fit the model
#' fit <- fit_binomial_glmm(y = y, X = X, Z = Z)
#' 
#' # 3. Print the summary
#' print(fit)
#' }
fit_binomial_glmm <- function(y, X, Z, initial_betas = NULL, initial_logsigma_u = 0) {
  
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
  cat("Fit from 'tmbayes' package\n\n") 
  cat("Formula: Binomial GLMM with random intercepts\n")
  cat("Optimizer status:", x$opt$message, "\n\n")
  
  cat("Parameter Estimates:\n")
  print(round(x$summary, 4))
  cat("\n")
}
