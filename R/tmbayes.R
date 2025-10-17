#' @title Fit a GLMM using a formula interface
#' @description This is the main user-facing function of the tmbayes package.
#' It allows fitting Generalized Linear Mixed Models (GLMMs) using a standard R formula.
#'
#' @param formula A two-sided linear formula object describing the fixed-effects and
#'   random-effects parts of the model. The structure is `response ~ fixed_effects + (random_effects | grouping_factor)`.
#' @param data A data frame containing the variables named in the formula.
#' @param family Currently, only `family = "binomial"` is supported.
#'
#' @return An object of class `tmbayes_fit`.
#' @import Matrix                # <--- AÑADE ESTA LÍNEA
#' @importFrom stats model.response
#' @importFrom methods as
#' @export
#' @examples
#' \dontrun{
#' # Simulate data for the formula interface
#' set.seed(1234)
#' n_groups <- 10
#' n_per_group <- 30
#'
#' # Create a data frame, which is required for formula interfaces
#' sim_data <- data.frame(
#'   group = factor(rep(1:n_groups, each = n_per_group)),
#'   x = rnorm(n_groups * n_per_group)
#' )
#'
#' # Simulate random effects for intercepts and slopes
#' u_intercepts <- rnorm(n_groups, 0, 1.0)
#' u_slopes <- rnorm(n_groups, 0, 0.5)
#'
#' # Calculate the linear predictor to generate the response variable
#' # Fixed effects: intercept=-1.0, slope=1.5
#' eta <- (-1.0 + u_intercepts[sim_data$group]) +
#'        (1.5 + u_slopes[sim_data$group]) * sim_data$x
#'
#' p <- 1 / (1 + exp(-eta))
#' sim_data$y <- rbinom(nrow(sim_data), size = 1, prob = p)
#'
#' # --- Fit the model using the new, simple formula interface! ---
#' fit <- tmbayes(
#'   formula = y ~ x + (x | group),
#'   data = sim_data,
#'   family = "binomial"
#' )
#'
#' print(fit)
#' }
tmbayes <- function(formula, data, family = "binomial") {

  # --- 1. Initial Validations ---
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("The 'formula' argument must be a valid R formula.")
  }
  if (missing(data)) {
    stop("The 'data' argument (a data frame) must be provided.")
  }
  if (family != "binomial") {
    stop("Currently, only the 'binomial' family is supported.")
  }
  # Check if lme4 is installed, as it's crucial
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("The 'lme4' package is required. Please install it with: install.packages('lme4')")
  }

  # --- 2. Formula Parsing with lme4 ---
  # lme4::lFormula does all the heavy lifting of parsing the mixed-effects formula.
  lmod <- lme4::lFormula(formula = formula, data = data)

  # Extract the components we need
  y <- stats::model.response(lmod$fr) # Response vector
  X <- lmod$X                        # Fixed-effects matrix
  
  # lme4's Z matrix is a transposed sparse matrix (Zt).
  # We convert it to a standard dense matrix to pass to TMB.
  Z <- as(t(lmod$reTrms$Zt), "matrix")

  # --- 3. Call the model fitting engine ---
  # Now that we have y, X, and Z, we just call your original function.
  fit <- fit_binomial_glmm(y = y, X = X, Z = Z)

  # --- 4. Return the result ---
  # The output object already has the correct class and print method.
  return(fit)
}
