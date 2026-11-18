#' Fit a latent cluster mixture model across strata
#'
#' Fits a finite overfitted latent cluster mixture model for a multivariate
#' covariate distribution observed in multiple strata. Continuous covariates
#' are standardized by default before calling the C++ sampler so that the
#' default component priors have a stable interpretation across variables.
#'
#' @param W_cont Matrix or data frame of continuous covariates. At least one
#'   continuous covariate is currently required by the C++ backend.
#' @param strata Vector giving the stratum membership of each observation.
#'   Arbitrary labels are accepted and internally recoded to 1,...,M.
#' @param W_cat Optional matrix or data frame of categorical covariates. Each
#'   column is internally converted to a factor and encoded as 1,...,C.
#' @param n_components Number of mixture components in the overfitted finite
#'   mixture.
#' @param alpha Total Dirichlet concentration parameter. The C++ sampler uses
#'   alpha / n_components for each component.
#' @param standardize Logical; standardize continuous covariates before fitting.
#' @param mu_sd Prior standard deviation of component means on the fitted scale.
#' @param prec_shape,prec_rate Shape and rate of the Gamma prior for component
#'   precisions on the fitted scale.
#' @param initial_Z Optional initial component labels. If NULL, k-means is used
#'   on a simple initialization representation of the covariates.
#' @param nstart Number of k-means random starts when initial_Z is NULL.
#' @param nburn Number of burn-in MCMC iterations.
#' @param thin Thinning interval.
#' @param n_draws Number of retained posterior draws.
#' @param seed Optional random seed. The C++ sampler uses R's RNG, so the same
#'   seed controls both initialization and posterior sampling.
#' @param verbose Logical; print C++ sampler progress.
#'
#' @return An object of class \code{strataLCM_fit}. Posterior arrays returned by
#'   the C++ sampler are available directly as components such as \code{PP},
#'   \code{Mu}, \code{Sig}, \code{Counts}, and \code{CountsV}.
#' @export
lcm_fit <- function(W_cont,
                    strata,
                    W_cat = NULL,
                    n_components = 20L,
                    alpha = 0.1,
                    standardize = TRUE,
                    mu_sd = 2.5,
                    prec_shape = 3,
                    prec_rate = 2,
                    initial_Z = NULL,
                    nstart = 10L,
                    nburn = 4000L,
                    thin = 10L,
                    n_draws = 4000L,
                    seed = NULL,
                    verbose = FALSE) {

  call <- match.call()

  Wc <- .lcm_continuous_matrix(W_cont)
  n <- nrow(Wc)

  if (length(strata) != n) {
    stop("'strata' must have the same number of observations as W_cont.",
         call. = FALSE)
  }
  if (anyNA(strata)) {
    stop("'strata' contains missing values.", call. = FALSE)
  }

  strata_factor <- factor(strata, levels = unique(strata))
  strata_code <- as.integer(strata_factor)
  strata_levels <- levels(strata_factor)

  n_components <- .lcm_scalar_integer(n_components, "n_components", lower = 1L)
  if (n_components > n) {
    stop("'n_components' cannot exceed the number of observations.", call. = FALSE)
  }
  nstart <- .lcm_scalar_integer(nstart, "nstart", lower = 1L)
  nburn <- .lcm_scalar_integer(nburn, "nburn", lower = 0L)
  thin <- .lcm_scalar_integer(thin, "thin", lower = 1L)
  n_draws <- .lcm_scalar_integer(n_draws, "n_draws", lower = 1L)

  .lcm_positive_scalar(alpha, "alpha")
  .lcm_positive_scalar(mu_sd, "mu_sd")
  .lcm_positive_scalar(prec_shape, "prec_shape")
  .lcm_positive_scalar(prec_rate, "prec_rate")

  prep_cont <- .lcm_standardize_continuous(Wc, standardize = standardize)
  Wc_fit <- prep_cont$x

  prep_cat <- .lcm_encode_categorical(W_cat, n_expected = n)

  if (!is.null(seed)) {
    seed <- .lcm_scalar_integer(seed, "seed", lower = 0L)
    set.seed(seed)
  }

  if (is.null(initial_Z)) {
    X_init <- .lcm_initialization_matrix(Wc_fit, prep_cat)
    initial_Z <- .lcm_kmeans_initialization(
      X_init,
      n_components = n_components,
      nstart = nstart
    )
  } else {
    if (length(initial_Z) != n) {
      stop("'initial_Z' must have the same number of observations as W_cont.",
           call. = FALSE)
    }
    if (anyNA(initial_Z) || any(!is.finite(initial_Z)) ||
        any(abs(initial_Z - round(initial_Z)) > 1e-8) ||
        any(initial_Z < 1) || any(initial_Z > n_components)) {
      stop("'initial_Z' must contain integer labels in 1,...,n_components.",
           call. = FALSE)
    }
    initial_Z <- as.integer(round(initial_Z))
  }

  if (is.null(prep_cat$x)) {
    ans <- LCM_cpp1(
      W_cont = Wc_fit,
      v = strata_code,
      initial_Z = initial_Z,
      mK = n_components,
      Alpha = alpha,
      mu_sd = mu_sd,
      prec_shape = prec_shape,
      prec_rate = prec_rate,
      NBURN = nburn,
      NTHIN = thin,
      NOUTSAMPLE = n_draws,
      print = as.integer(isTRUE(verbose))
    )
  } else {
    ans <- LCM_cpp(
      W_cont = Wc_fit,
      W_cat = prep_cat$x,
      v = strata_code,
      initial_Z = initial_Z,
      mK = n_components,
      Alpha = alpha,
      mu_sd = mu_sd,
      prec_shape = prec_shape,
      prec_rate = prec_rate,
      NBURN = nburn,
      NTHIN = thin,
      NOUTSAMPLE = n_draws,
      print = as.integer(isTRUE(verbose))
    )
  }

  ans$call <- call
  ans$strata_code <- strata_code
  ans$strata_levels <- strata_levels
  ans$stratum_sizes <- tabulate(strata_code, nbins = length(strata_levels))
  ans$initial_Z <- initial_Z
  ans$continuous <- list(
    names = colnames(Wc),
    standardized = isTRUE(standardize),
    center = prep_cont$center,
    scale = prep_cont$scale
  )
  ans$categorical <- list(
    names = prep_cat$names,
    levels = prep_cat$levels
  )
  ans$settings <- list(
    n_components = n_components,
    alpha = alpha,
    mu_sd = mu_sd,
    prec_shape = prec_shape,
    prec_rate = prec_rate,
    nburn = nburn,
    thin = thin,
    n_draws = n_draws,
    seed = seed
  )

  class(ans) <- c("strataLCM_fit", "list")
  ans
}

#' @export
print.strataLCM_fit <- function(x, ...) {
  cat("strataLCM fit\n")
  cat("  observations :", sum(x$stratum_sizes), "\n")
  cat("  strata       :", length(x$stratum_sizes), "\n")
  cat("  components   :", x$settings$n_components, "\n")
  cat("  draws        :", x$settings$n_draws, "\n")
  cat("  continuous   :", length(x$continuous$names), "\n")
  cat("  categorical  :", length(x$categorical$names), "\n")
  cat("  standardized :", x$continuous$standardized, "\n")
  invisible(x)
}

.lcm_continuous_matrix <- function(x) {
  if (is.null(x)) {
    stop("W_cont is required by the current C++ backend.", call. = FALSE)
  }
  x <- as.data.frame(x, check.names = FALSE)
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop("W_cont must have at least one row and one column.", call. = FALSE)
  }
  ok <- vapply(x, is.numeric, logical(1))
  if (!all(ok)) {
    stop("All columns of W_cont must be numeric.", call. = FALSE)
  }
  out <- as.matrix(x)
  storage.mode(out) <- "double"
  if (anyNA(out) || any(!is.finite(out))) {
    stop("W_cont contains NA/NaN/Inf.", call. = FALSE)
  }
  if (is.null(colnames(out))) {
    colnames(out) <- paste0("W", seq_len(ncol(out)))
  }
  out
}

.lcm_standardize_continuous <- function(x, standardize) {
  p <- ncol(x)
  if (!isTRUE(standardize)) {
    return(list(
      x = x,
      center = rep(0, p),
      scale = rep(1, p)
    ))
  }

  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  bad <- !is.finite(scale) | scale <= 0
  if (any(bad)) {
    stop(
      "Cannot standardize constant/non-varying continuous covariate(s): ",
      paste(colnames(x)[bad], collapse = ", "),
      call. = FALSE
    )
  }

  xs <- sweep(x, 2L, center, FUN = "-")
  xs <- sweep(xs, 2L, scale, FUN = "/")
  colnames(xs) <- colnames(x)
  list(x = xs, center = center, scale = scale)
}

.lcm_encode_categorical <- function(x, n_expected) {
  if (is.null(x)) {
    return(list(x = NULL, names = character(0), levels = list()))
  }

  df <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) != n_expected) {
    stop("W_cat must have the same number of observations as W_cont.",
         call. = FALSE)
  }
  if (ncol(df) < 1L) {
    return(list(x = NULL, names = character(0), levels = list()))
  }
  if (anyNA(df)) {
    stop("W_cat contains missing values.", call. = FALSE)
  }
  if (is.null(names(df))) {
    names(df) <- paste0("C", seq_len(ncol(df)))
  }

  levs <- vector("list", ncol(df))
  out <- matrix(NA_real_, nrow(df), ncol(df))
  colnames(out) <- names(df)

  for (j in seq_len(ncol(df))) {
    z <- df[[j]]
    if (is.factor(z)) {
      f <- droplevels(z)
    } else {
      # factor() gives a deterministic mapping; the original levels are stored.
      f <- factor(z)
    }
    if (nlevels(f) < 2L) {
      stop("Categorical covariate '", names(df)[j],
           "' has fewer than two observed levels.", call. = FALSE)
    }
    levs[[j]] <- levels(f)
    out[, j] <- as.integer(f)
  }

  list(x = out, names = names(df), levels = levs)
}

.lcm_initialization_matrix <- function(W_cont_fit, prep_cat) {
  X <- W_cont_fit
  if (!is.null(prep_cat$x)) {
    cat_parts <- vector("list", ncol(prep_cat$x))
    for (j in seq_len(ncol(prep_cat$x))) {
      Cj <- length(prep_cat$levels[[j]])
      tmp <- matrix(0, nrow(prep_cat$x), Cj)
      tmp[cbind(seq_len(nrow(tmp)), prep_cat$x[, j])] <- 1
      colnames(tmp) <- paste0(prep_cat$names[j], "_", seq_len(Cj))
      cat_parts[[j]] <- tmp
    }
    X <- cbind(X, do.call(cbind, cat_parts))
  }
  X
}

.lcm_kmeans_initialization <- function(X, n_components, nstart) {
  n_unique <- nrow(unique(as.data.frame(X)))
  centers <- min(n_components, n_unique)
  if (centers <= 1L) {
    return(rep.int(1L, nrow(X)))
  }
  out <- stats::kmeans(
    X,
    centers = centers,
    nstart = nstart,
    iter.max = 100L
  )$cluster
  as.integer(out)
}

.lcm_positive_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop("'", name, "' must be one finite positive number.", call. = FALSE)
  }
  invisible(x)
}

.lcm_scalar_integer <- function(x, name, lower = 0L) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) ||
      abs(x - round(x)) > 1e-8 || x < lower) {
    stop("'", name, "' must be one integer >= ", lower, ".", call. = FALSE)
  }
  as.integer(round(x))
}
