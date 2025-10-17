# tmbayes <img src="httpslogo.png" align="right" height="120" /> 

`tmbayes` is an R package designed to make the powerful **Template Model Builder (TMB)** engine accessible and easy to use for applied statisticians and researchers. It provides user-friendly wrappers for common statistical models, allowing you to leverage TMB's speed and flexibility without writing or compiling C++ code yourself.

## The Motivation

TMB is an outstanding tool for fitting complex statistical models, offering performance comparable to established platforms like ADMB. However, its reliance on C++ templates can present a steep learning curve for those who primarily work in R. The process of writing the model syntax in C++, compiling it, and linking it within an R session can be a significant barrier.

`tmbayes` bridges this gap by providing pre-compiled TMB models that can be called directly from R with simple, intuitive functions.

## Core Idea: Laplace Approximation for Bayesian-like Inference

This package uses the **Laplace Approximation**, a method that provides a fast and accurate Gaussian approximation to the posterior distribution of model parameters. This gives results that are "Bayesian-like" in their interpretation but are obtained with the speed of likelihood optimization. The name `tmbayes` reflects this connection between TMB and the Bayesian-esque nature of the Laplace Approximation.

## Installation

You can install the development version of `tmbayes` from [GitHub](https://github.com/) with:

```r
# First, ensure you have the devtools package
# install.packages("devtools")

# Now, install tmbayes
devtools::install_github("giovannitasso/tmbayes")
```

**Important Note:** Since this package compiles C++ code, you will need the appropriate development tools for your operating system:
* **Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
* **macOS:** Install the Command Line Tools by running `xcode-select --install` in your terminal.
* **Linux (Debian/Ubuntu):** Install `r-base-dev` and `build-essential` by running `sudo apt-get install r-base-dev build-essential`.

## Quick Start: Fitting a Binomial GLMM

Here is a complete example of how to simulate data and fit a Generalized Linear Mixed Model (GLMM) for a binary outcome.

```r
# 1. Load the library
library(tmbayes)

# 2. Simulate data
# We'll simulate data where the effect of predictor 'x' varies across groups.
set.seed(1234)
n_groups <- 10
n_per_group <- 30 # More data helps estimate random slopes
n_tot <- n_groups * n_per_group

group <- rep(1:n_groups, each = n_per_group)
x <- rnorm(n_tot) 

# True fixed effects
beta0_true <- -1.0  # Intercept
beta1_true <- 1.5   # Average slope for x

# True random effects standard deviations
sigma_intercept_true <- 1.0 # Variability of intercepts across groups
sigma_slope_true <- 0.5     # Variability of slopes across groups

# Simulate random effects for each group from a Normal distribution
u_intercepts <- rnorm(n_groups, 0, sigma_intercept_true) # Random intercepts
u_slopes <- rnorm(n_groups, 0, sigma_slope_true)       # Random slopes

# Calculate linear predictor with group-specific intercepts and slopes
eta <- (beta0_true + u_intercepts[group]) + (beta1_true + u_slopes[group]) * x
p <- 1 / (1 + exp(-eta))
y <- rbinom(n_tot, size = 1, prob = p)

# 3. Define design matrices
# Fixed effects matrix (X)
X <- cbind(1, x)
colnames(X) <- c("Intercept", "x_variable")

# Random effects matrix (Z) for intercepts and slopes
# This requires manual construction. The vector of random effects 'u' in TMB will be
# structured as [intercept_1, ..., intercept_10, slope_1, ..., slope_10].
Z <- matrix(0, nrow = n_tot, ncol = n_groups * 2)

for (i in 1:n_tot) {
  group_idx <- group[i]
  # Part 1: Assign the random intercept effect
  Z[i, group_idx] <- 1 
  # Part 2: Assign the random slope effect, multiplied by the predictor value 'x'
  Z[i, n_groups + group_idx] <- x[i]
}

# 4. Fit the model
# Our C++ model assumes all 20 random effects (10 intercepts + 10 slopes)
# are drawn from a single distribution with one standard deviation, sigma_u.
# This is an approximation of the true data-generating process.
fit <- fit_binomial_glmm(y = y, X = X, Z = Z)

# 5. Print the results
print(fit)
```

### Expected Output

```
Laplace Approximation Fit from 'tmbayes'

Formula: Binomial GLMM with random intercepts
Optimizer status: relative convergence (4) 

Parameter Estimates:
             Estimate Std. Error  z value Pr(>|z|)
Intercept     -1.0711     0.3705  -2.8911   0.0038
x_variable     1.6521     0.2285   7.2307   0.0000
sigma_u        0.9829     0.2889   3.4020   0.0007
```

## Features

* **No C++ Required:** Fit complex TMB models using only familiar R functions.
* **Fast Estimation:** Leverages TMB's high-performance C++ backend for speed.
* **Simple Interface:** Designed to be intuitive for users familiar with R's standard modeling syntax (like `lm` or `lme4`).
* **Extensible:** Provides a clear framework for adding new pre-compiled TMB models in the future.

## Available Models

Currently, `tmbayes` includes:

* `fit_binomial_glmm()`: Fits a Generalized Linear Mixed Model for binary outcomes (0/1) with a random intercept for grouping factors.

## Contributing

Feedback, bug reports, and feature requests are welcome! Please open an issue on the [GitHub repository issues page](https://github.com/giovannitasso/tmbayes/issues).
