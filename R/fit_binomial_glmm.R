#' @title Fit a Binomial Generalized Linear Mixed Model with TMB (Simple Version)
#' @description This function fits a binomial GLMM using a pre-compiled TMB model.
#' This version uses a simplified random effects structure (single variance component).
#'
#' @param y A numeric vector representing the binary response variable (0s and 1s).
#' @param X A numeric matrix of fixed effects covariates. An intercept is recommended.
#' @param Z A numeric matrix for the random effects design.
#' @param initial_betas A numeric vector for the starting values of the fixed effects (betas).
#'        If NULL, defaults to zeros.
#' @param initial_logsigma_u A single numeric value for the starting value of the log standard
#'        deviation of the random effects. Defaults to 0.
#'
#' @return A list object of class `tmbr_fit` containing the model results.
#' @export
#' @examples
#' \dontrun{
#' # This function is typically called internally by tmbr().
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
    # n_groups removed for this simple version
  )
  
  # ---- 2. Parameter Initialization ----
  if (is.null(initial_betas)) {
    initial_betas <- rep(0, ncol(X))
  }
  
  # Reverted to simple parameterization
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
    DLL = "tmbrmodels",
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
  
  # 1. Fixed Effects
  num_fixed_effects <- ncol(X)
  fixed_eff_idx <- which(rownames(summary_report) == "betas")
  correct_fixed_eff_idx <- fixed_eff_idx[1:num_fixed_effects]
  fixed_coeffs <- summary_report[correct_fixed_eff_idx, , drop = FALSE]
  rownames(fixed_coeffs) <- if (!is.null(colnames(X))) colnames(X) else paste0("beta_", 1:nrow(fixed_coeffs))
  
  # 2. Random Effects Variance Component (Simple Version)
  random_eff_idx <- which(rownames(summary_report) == "sigma_u")
  random_coeffs <- summary_report[random_eff_idx, , drop = FALSE]
  rownames(random_coeffs) <- "sigma_u"

  # 3. Combine and Format
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
  
  class(output) <- "tmbr_fit" 
  
  return(output)
}


#' @title Print method for tmbr_fit objects
#' @description A custom print method to display a concise summary of the fitted model.
#' @param x An object of class `tmbr_fit`.
#' @param ... Additional arguments passed to `print`.
#' @method print tmbr_fit
#' @export
print.tmbr_fit <- function(x, ...) {

  cat("Laplace Approximation Fit from 'tmbrmodels'\n\n") 
  

  cat("Formula: Binomial GLMM\n")
  cat("Optimizer status:", x$opt$message, "\n\n")
  
  cat("Parameter Estimates:\n")
  print(round(x$summary, 4))
  cat("\n")
}
