# strataLCM 0.1.0

Initial public release.

* Fit finite latent cluster mixture models with shared components and
  stratum-specific mixing weights using `lcm_fit()`.
* Support continuous covariates with Gaussian kernels and categorical
  covariates with component-specific categorical distributions.
* Generate posterior-predictive covariates for selected strata and MCMC draws
  using `lcm_draw_covariates()`.
* Summarize component-based support across strata using `lcm_support()`,
  including external support and mass assigned to globally empty components.
* Provide reproducible R-based random-number generation in the MCMC sampler.
