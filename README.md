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

# 2. Simulate data for a mixed-effects model
# We'll simulate data from 10 groups, with 20 observations per group.
set.seed(1234)
n_groups <- 10
n_per_group <- 20

group <- rep(1:n_groups, each = n_per_group) # Grouping variable
x <- rnorm(n_groups * n_per_group)          # A fixed-effect predictor

# True parameter values
beta0_true <- -1.0  # Intercept
beta1_true <- 1.5   # Slope for x
sigma_u_true <- 1.0 # Standard deviation of random effects

# Simulate random intercepts for each group
u_true <- rnorm(n_groups, 0, sigma_u_true)

# Calculate the linear predictor and probabilities
eta <- beta0_true + beta1_true * x + u_true[group]
p <- 1 / (1 + exp(-eta))
y <- rbinom(n_groups * n_per_group, size = 1, prob = p) # Binary response

# 3. Define design matrices
# The 'X' matrix is for fixed effects (Intercept and x)
X <- cbind(1, x)
colnames(X) <- c("Intercept", "x_variable")

# The 'Z' matrix is for random effects (a random intercept for each group)
Z <- model.matrix(~ factor(group) - 1)

# 4. Fit the model using the package's main function
# All the C++ compilation and TMB machinery is handled internally here.
fit <- fit_binomial_glmm(y = y, X = X, Z = Z)

# 5. Print the results
# The custom print method provides a clean and readable summary.
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
