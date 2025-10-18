<p align="center">
  <img src="man/figures/logotmbayes.png" height="130" />
</p>

# tmbayes

`tmbayes` is an R package designed to make the powerful **Template Model Builder (TMB)** engine accessible and easy to use for applied statisticians and researchers. It provides a high-level formula interface for common statistical models, allowing you to leverage TMB's speed and flexibility without writing or compiling C++ code yourself.


## The Motivation

TMB is an outstanding tool for fitting complex statistical models. However, its reliance on C++ templates can present a steep learning curve for those who primarily work in R. The process of writing the model syntax in C++, compiling it, and linking it within an R session can be a significant barrier.

`tmbayes` bridges this gap by providing a simple, formula-based interface that handles all the C++ and TMB complexity behind the scenes.


## Core Idea: Laplace Approximation for Bayesian-like Inference

This package uses the **Laplace Approximation**, a method that provides a fast and accurate Gaussian approximation to the posterior distribution of model parameters. This gives results that are "Bayesian-like" in their interpretation but are obtained with the speed of likelihood optimization. The name `tmbayes` reflects this connection between TMB and the Bayesian-esque nature of the Laplace Approximation.


## Installation

You can install the development version of `tmbayes` from GitHub with:

```r
# First, ensure you have the devtools package
# install.packages("devtools")

# Now, install tmbayes
# This will also install its dependencies, like TMB and lme4.
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
# We'll create a data frame, which is the standard for formula-based modeling in R.
set.seed(1234)
n_groups <- 10
n_per_group <- 30

sim_data <- data.frame(
  group = factor(rep(1:n_groups, each = n_per_group)),
  x = rnorm(n_groups * n_per_group)
)

# True parameter values
fixed_intercept <- -1.0
fixed_slope <- 1.5
sigma_intercept <- 1.0 # SD of random intercepts
sigma_slope <- 0.5     # SD of random slopes

# Simulate random effects for each group
random_intercepts <- rnorm(n_groups, 0, sigma_intercept)
random_slopes <- rnorm(n_groups, 0, sigma_slope)

# Calculate the linear predictor and response
eta <- (fixed_intercept + random_intercepts[sim_data$group]) + 
       (fixed_slope + random_slopes[sim_data$group]) * sim_data$x
       
p <- 1 / (1 + exp(-eta))
sim_data$y <- rbinom(nrow(sim_data), size = 1, prob = p)

# 3. Fit the model using the simple formula interface
# The complexity of creating design matrices is handled automatically!
fit <- tmbayes(
  formula = y ~ x + (x | group), 
  data = sim_data,
  family = "binomial"
)

# 4. Print the results
print(fit)

```

### Expected Output

```
Laplace Approximation Fit from 'tmbayes'

Formula: Binomial GLMM
Optimizer status: relative convergence (4) 

Parameter Estimates:
             Estimate Std. Error  z value Pr(>|z|)
(Intercept)   -1.1259     0.3649  -3.0855   0.0020
x              1.6083     0.2224   7.2311   0.0000
sigma_u        0.8915     0.1652   5.3965   0.0000
```
*Note: The estimated `sigma_u` represents the standard deviation for all random effects combined, as the current model uses a single variance component.*


## Features

* **Formula Interface:** Fit complex models using the standard R formula syntax (e.g., `y ~ x + (x | group)`), just like in `lme4`.
* **No C++ Required:** All the TMB C++ code is pre-built and handled internally.
* **Fast Estimation:** Leverages TMB's high-performance C++ backend for speed.
* **Extensible:** Provides a clear framework for adding new model families in the future.


## Available Models

Currently, `tmbayes` provides the main `tmbayes()` function, which supports:

* **Binomial GLMMs:** For binary outcomes (0/1).
* **Flexible Random Effects:** Supports random intercepts `(1 | group)`, random slopes `(0 + x | group)`, and both combined `(x | group)`.


## Contributing

Feedback, bug reports, and feature requests are welcome! Please open an issue on the [GitHub repository issues page](https://github.com/giovannitasso/tmbayes/issues).


