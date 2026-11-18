# strataLCM

**Latent Cluster Mixture Models for Stratum-Specific Covariate Distributions**

`strataLCM` is an R package for Bayesian estimation of multivariate covariate distributions across predefined strata. It uses a finite latent cluster mixture (LCM) model with **shared mixture components** and **stratum-specific mixture weights**, allowing information to be shared across strata while retaining differences in their covariate distributions.

The package supports continuous and categorical covariates, posterior-predictive covariate generation, and diagnostics describing how shared mixture components connect strata. It is designed to provide estimated covariate distributions that can be used in downstream analyses, including outcome-model standardization.

## Installation

The package is not currently distributed through CRAN. After the GitHub repository is published, it can be installed with:

```r
install.packages("remotes")
remotes::install_github("tom-ohigashi/strataLCM")
```

The package uses `Rcpp` and `RcppArmadillo` and requires a working C++ build toolchain for installation from source.

## Model overview

For stratum \(v\), the covariate distribution is modeled as

$$
f_v(w) = \sum_{\ell=1}^{L} \pi_{v\ell}\,h_\ell(w),
$$

where \(h_\ell\) denotes a mixture component shared across strata and \(\pi_{v\ell}\) is its stratum-specific weight. Within each component, covariates are modeled as conditionally independent: continuous variables have Gaussian kernels, and categorical variables have categorical probability distributions.

The mixture weights have independent symmetric Dirichlet priors,

$$
\boldsymbol{\pi}_v \sim
\operatorname{Dirichlet}(\alpha/L,\ldots,\alpha/L),
$$

where `alpha` is the **total concentration parameter**, not the concentration parameter of each individual component. The model is fitted by MCMC using an `RcppArmadillo` computational backend.

## Quick start

The following example fits an LCM to three strata using two continuous covariates and two categorical covariates.

```r
library(strataLCM)

set.seed(123)

strata <- rep(c("A", "B", "C"), times = c(80, 40, 15))
n <- length(strata)

W_cont <- data.frame(
  age = rnorm(n, mean = c(A = 50, B = 57, C = 62)[strata], sd = 9),
  biomarker = rnorm(n, mean = c(A = 0, B = 1, C = -1)[strata], sd = 1)
)

W_cat <- data.frame(
  smoker = sample(c("no", "yes"), n, replace = TRUE),
  site = sample(c("X", "Y", "Z"), n, replace = TRUE)
)

fit <- lcm_fit(
  W_cont = W_cont,
  W_cat = W_cat,
  strata = strata,
  n_components = 12,
  alpha = 0.1,
  nburn = 200,
  thin = 2,
  n_draws = 100,
  seed = 2026
)

print(fit)
```

The relatively short MCMC run above is for demonstration only. For substantive analyses, choose the number of components and MCMC settings carefully and assess posterior stability and convergence.

### Generate posterior-predictive covariates

Generate 200 covariate vectors for stratum `"C"`, conditional on the fifth retained posterior draw:

```r
W_new <- lcm_draw_covariates(
  fit,
  stratum = "C",
  n = 200,
  draws = 5,
  seed = 2027
)

head(W_new)
```

The function returns continuous covariates on their original scale by default, and categorical covariates as factors with their original levels. With a single value of `draws`, the result is a data frame; with multiple posterior draws, it is a named list of data frames.

These generated covariates can be passed to a separately fitted outcome regression model for standardization. **Outcome regression fitting and treatment-effect estimation are not currently implemented as integrated package functions.**

### Examine information-sharing diagnostics

```r
support <- lcm_support(fit)

round(support$mean, 3)
support$external_support_mean
support$empty_mass_mean
```

`support$mean` is a matrix whose rows correspond to target strata and columns to strata contributing observations to shared mixture components. `external_support_mean` summarizes off-diagonal support, while `empty_mass_mean` reports posterior mass assigned to components with no allocated observations.

By default, support is **not renormalized** after excluding globally empty components, so rows of `support$mean` need not sum to one. Set `conditional = TRUE` to normalize support over occupied components within each posterior draw:

```r
support_conditional <- lcm_support(fit, conditional = TRUE)
round(support_conditional$mean, 3)
```

These are model-based support diagnostics, **not effective sample sizes**, and should not be interpreted as direct evidence of causal identification or improved estimation accuracy.

## Main functions

| Function | Purpose |
| --- | --- |
| `lcm_fit()` | Fit a finite LCM with stratum-specific mixture weights and shared component parameters. |
| `lcm_draw_covariates()` | Generate posterior-predictive covariates for a specified stratum and selected MCMC draws. |
| `lcm_support()` | Summarize component-based support attributed to observed strata. |

Lower-level C++ sampler interfaces, `LCM_cpp()` and `LCM_cpp1()`, are also exported. Most users should prefer `lcm_fit()`, which handles stratum labels, categorical coding, and optional standardization.

## Input conventions and current limitations

- `W_cont` must contain at least one numeric continuous covariate; models containing only categorical covariates are not currently supported by the fitting backend.
- `W_cat` is optional. Each categorical variable should be a separate column; `lcm_fit()` internally encodes observed levels, including 0/1 and character-valued categories.
- Missing or non-finite continuous values and missing categorical values are not accepted. Handle missingness before fitting.
- Continuous covariates are standardized by default (`standardize = TRUE`). Posterior-predictive draws are automatically transformed back to the original scale unless `original_scale = FALSE` is requested.
- Strata must be defined in advance. The model uses a fixed, user-specified number of mixture components and assumes conditional independence of covariates within components.
- Numerical stability, sensitivity to priors and the number of components, and MCMC convergence should be examined for each application.

For detailed arguments and returned objects, see `?lcm_fit`, `?lcm_draw_covariates`, and `?lcm_support` after installation.

## Development status

This repository provides the initial `0.1.0` release of `strataLCM`. The current focus is estimation of stratum-specific covariate distributions and their posterior-predictive use; the API and documentation may evolve.

## License

GPL-3.
