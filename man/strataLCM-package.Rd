
\name{strataLCM-package}
\alias{strataLCM-package}
\alias{strataLCM}
\docType{package}

\title{Latent Cluster Mixture Models for Stratum-Specific Covariates}

\description{
The \code{strataLCM} package provides Bayesian finite mixture
models for estimating multivariate covariate distributions
within predefined strata.
}

\details{
The model uses shared mixture components and stratum-specific
mixing weights to allow information borrowing across strata
while accommodating distributional heterogeneity.

The package supports continuous and categorical covariates,
posterior-predictive covariate generation, and diagnostics
of information sharing through shared mixture components.
}

\seealso{
\code{\link{lcm_fit}},
\code{\link{lcm_draw_covariates}},
\code{\link{lcm_support}}
}

\keyword{package}
