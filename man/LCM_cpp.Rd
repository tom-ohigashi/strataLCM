\name{LCM_cpp}
\alias{LCM_cpp}
\alias{LCM_cpp1}
\title{Latent cluster mixture MCMC samplers}
\description{
Low-level RcppArmadillo samplers for stratum-specific latent cluster mixture models.
\code{LCM_cpp} handles continuous and categorical covariates; \code{LCM_cpp1}
handles continuous covariates only.
}
\usage{
LCM_cpp(W_cont, W_cat, v, initial_Z, mK = 7L, Alpha = 0.1,
  mu_sd = 1000, prec_shape = 0.1, prec_rate = 0.1,
  NBURN = 4000L, NTHIN = 10L,
  NOUTSAMPLE = 4000L, print = 0L)

LCM_cpp1(W_cont, v, initial_Z, mK = 7L, Alpha = 0.1,
  mu_sd = 1000, prec_shape = 0.1, prec_rate = 0.1,
  NBURN = 4000L, NTHIN = 10L,
  NOUTSAMPLE = 4000L, print = 0L)
}
\arguments{
  \item{W_cont}{Numeric matrix of continuous covariates.}
  \item{W_cat}{Numeric matrix of categorical covariates coded as consecutive integers 1, ..., C for each column.}
  \item{v}{Stratum labels coded as consecutive integers 1, ..., M.}
  \item{initial_Z}{Initial mixture allocation labels coded 1, ..., mK.}
  \item{mK}{Maximum number of mixture components.}
  \item{Alpha}{Total Dirichlet concentration for stratum-specific mixture weights.}
  \item{mu_sd}{Prior standard deviation for component means.}
  \item{prec_shape}{Gamma prior shape for component precisions.}
  \item{prec_rate}{Gamma prior rate for component precisions.}
  \item{NBURN}{Number of burn-in iterations.}
  \item{NTHIN}{Thinning interval.}
  \item{NOUTSAMPLE}{Number of retained posterior draws.}
  \item{print}{Set to 1 for progress messages.}
}
\value{
A list containing posterior arrays \code{Mu}, \code{Sig}, \code{PP}, and
\code{Counts}; \code{LCM_cpp} additionally returns categorical probabilities
\code{Pi}. Both functions also return a \code{Hyperparameters} list.
}
\details{
This is a low-level prototype interface. In the mixed-data sampler, categorical
variables must not be dummy-expanded and 0/1 variables must be recoded as 1/2.
All random-number generation in the revised sampler uses R's RNG, so results are
reproducible under \code{set.seed()}.
}
