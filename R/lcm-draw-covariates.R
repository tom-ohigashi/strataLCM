#' Draw covariates from a fitted stratum-specific LCM
#'
#' Generates posterior-predictive covariate draws from a selected stratum and
#' one or more retained MCMC draws. The component label is first sampled from
#' the stratum-specific mixing weights, after which continuous and categorical
#' covariates are sampled from the selected shared component.
#'
#' Continuous covariates are returned on their original scale by default. If
#' the model was fitted after standardization, the stored centers and scales
#' are used to transform posterior-predictive draws back automatically.
#'
#' @param object A fitted \code{strataLCM_fit} object.
#' @param stratum Stratum label for which posterior-predictive covariates are
#'   required. Use the same label as supplied to \code{lcm_fit()}.
#' @param n Number of posterior-predictive observations to generate for each
#'   selected posterior draw.
#' @param draws Integer vector indexing retained MCMC draws. By default all
#'   retained draws are used.
#' @param original_scale Logical. If TRUE (default), continuous covariates are
#'   returned on the original data scale; if FALSE, they are returned on the
#'   fitted (usually standardized) scale.
#' @param include_component Logical. If TRUE, add the sampled component label
#'   as column \code{.component}.
#' @param seed Optional random seed.
#'
#' @return If one posterior draw is requested, a data frame with \code{n} rows.
#'   If multiple posterior draws are requested, a named list of data frames,
#'   one for each retained posterior draw. Each data frame has attributes
#'   \code{posterior_draw} and \code{stratum}.
#' @export
lcm_draw_covariates <- function(object,
                                stratum,
                                n = 1000L,
                                draws = NULL,
                                original_scale = TRUE,
                                include_component = FALSE,
                                seed = NULL) {
  if (!inherits(object, "strataLCM_fit")) {
    stop("'object' must inherit from 'strataLCM_fit'.", call. = FALSE)
  }

  n <- .lcm_scalar_integer(n, "n", lower = 1L)

  if (length(stratum) != 1L || is.na(stratum)) {
    stop("'stratum' must be one non-missing stratum label.", call. = FALSE)
  }

  stratum_chr <- as.character(stratum)
  v <- match(stratum_chr, object$strata_levels)
  if (is.na(v)) {
    stop(
      "Unknown stratum label. Available strata are: ",
      paste(object$strata_levels, collapse = ", "),
      call. = FALSE
    )
  }

  PP <- object$PP
  Mu <- object$Mu
  Sig <- object$Sig

  if (length(dim(PP)) != 3L || length(dim(Mu)) != 3L ||
      length(dim(Sig)) != 3L) {
    stop("Malformed fitted object: PP, Mu, and Sig must be 3-dimensional arrays.",
         call. = FALSE)
  }

  S <- dim(PP)[3L]
  L <- dim(PP)[2L]
  p_cont <- dim(Mu)[1L]

  if (!identical(dim(Mu), dim(Sig)) || dim(Mu)[2L] != L || dim(Mu)[3L] != S) {
    stop("Malformed fitted object: PP, Mu, and Sig dimensions are inconsistent.",
         call. = FALSE)
  }

  if (is.null(draws)) {
    draws <- seq_len(S)
  } else {
    if (anyNA(draws) || any(!is.finite(draws)) ||
        any(abs(draws - round(draws)) > 1e-8)) {
      stop("'draws' must contain finite integer indices.", call. = FALSE)
    }
    draws <- as.integer(round(draws))
    if (any(draws < 1L | draws > S)) {
      stop("'draws' contains an index outside 1,...,number of posterior draws.",
           call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    seed <- .lcm_scalar_integer(seed, "seed", lower = 0L)
    set.seed(seed)
  }

  cont_names <- object$continuous$names
  if (length(cont_names) != p_cont) {
    cont_names <- paste0("W", seq_len(p_cont))
  }

  center <- object$continuous$center
  scale <- object$continuous$scale
  if (length(center) != p_cont || length(scale) != p_cont) {
    stop("Malformed fitted object: stored continuous center/scale are inconsistent.",
         call. = FALSE)
  }

  p_cat <- length(object$categorical$levels)
  has_cat <- p_cat > 0L

  if (has_cat) {
    if (is.null(object$Pi) || length(object$Pi) != p_cat) {
      stop("Malformed fitted object: categorical posterior probabilities are missing.",
           call. = FALSE)
    }
  }

  draw_one <- function(m) {
    prob <- as.numeric(PP[v, , m])
    if (any(!is.finite(prob)) || any(prob < 0) || sum(prob) <= 0) {
      stop("Invalid component probabilities in posterior draw ", m, ".", call. = FALSE)
    }
    prob <- prob / sum(prob)

    z <- sample.int(L, size = n, replace = TRUE, prob = prob)

    Wc <- matrix(NA_real_, nrow = n, ncol = p_cont)
    colnames(Wc) <- cont_names

    # Generate by component to avoid an R loop over observations.
    for (ell in unique(z)) {
      ii <- which(z == ell)
      nell <- length(ii)
      for (j in seq_len(p_cont)) {
        Wc[ii, j] <- stats::rnorm(
          nell,
          mean = Mu[j, ell, m],
          sd = Sig[j, ell, m]
        )
      }
    }

    if (isTRUE(original_scale)) {
      Wc <- sweep(Wc, 2L, scale, FUN = "*")
      Wc <- sweep(Wc, 2L, center, FUN = "+")
      colnames(Wc) <- cont_names
    }

    out <- as.data.frame(Wc, check.names = FALSE)

    if (has_cat) {
      cat_names <- object$categorical$names
      if (length(cat_names) != p_cat) {
        cat_names <- paste0("C", seq_len(p_cat))
      }

      for (j in seq_len(p_cat)) {
        lev <- object$categorical$levels[[j]]
        Cj <- length(lev)
        xj <- integer(n)

        Pij <- object$Pi[[j]]
        if (length(dim(Pij)) != 3L || dim(Pij)[2L] != L || dim(Pij)[3L] != S ||
            dim(Pij)[1L] < Cj) {
          stop("Malformed fitted object: categorical probability array ", j,
               " has inconsistent dimensions.", call. = FALSE)
        }

        for (ell in unique(z)) {
          ii <- which(z == ell)
          pj <- as.numeric(Pij[seq_len(Cj), ell, m])
          if (any(!is.finite(pj)) || any(pj < 0) || sum(pj) <= 0) {
            stop("Invalid categorical probabilities in posterior draw ", m,
                 ", variable ", j, ", component ", ell, ".", call. = FALSE)
          }
          pj <- pj / sum(pj)
          xj[ii] <- sample.int(Cj, size = length(ii), replace = TRUE, prob = pj)
        }

        out[[cat_names[j]]] <- factor(lev[xj], levels = lev)
      }
    }

    if (isTRUE(include_component)) {
      out$.component <- z
    }

    attr(out, "posterior_draw") <- m
    attr(out, "stratum") <- object$strata_levels[v]
    out
  }

  ans <- lapply(draws, draw_one)
  names(ans) <- paste0("draw_", draws)

  if (length(ans) == 1L) ans[[1L]] else ans
}
