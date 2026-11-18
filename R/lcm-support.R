#' Posterior support/contribution matrices across strata
#'
#' Uses the saved stratum-by-component allocation counts to quantify which
#' observed strata support the mixture components used by each target stratum.
#'
#' For posterior draw m,
#' B[v,u,m] = sum_l pi[v,l,m] * n[u,l,m] / n[.,l,m]
#' over occupied components. Hence each row sums to one minus the posterior
#' mass assigned to globally empty components.
#'
#' @param object A fitted \code{strataLCM_fit} object containing CountsV.
#' @param conditional Logical. If TRUE, normalize each target-stratum row within
#'   each posterior draw over occupied components before summarizing.
#' @param keep_draws Logical. If TRUE, return the full K x K x S posterior array.
#'
#' @return A list containing the posterior mean support matrix, posterior mean
#'   empty-component mass, posterior mean external support, and optionally draws.
#' @export
lcm_support <- function(object, conditional = FALSE, keep_draws = FALSE) {
  if (is.null(object$CountsV)) {
    stop("The fitted object does not contain CountsV. Refit with the updated C++ sampler.",
         call. = FALSE)
  }

  PP <- object$PP
  CountsV <- object$CountsV

  K <- dim(CountsV)[1L]
  L <- dim(CountsV)[2L]
  S <- dim(CountsV)[3L]

  if (!identical(dim(PP), c(K, L, S))) {
    stop("PP and CountsV dimensions are inconsistent.", call. = FALSE)
  }

  B <- array(0, dim = c(K, K, S))
  empty_mass <- matrix(0, nrow = S, ncol = K)

  for (m in seq_len(S)) {
    n_l <- colSums(CountsV[, , m, drop = FALSE][, , 1L])
    occupied <- n_l > 0
    empty <- !occupied

    for (v in seq_len(K)) {
      if (any(empty)) {
        empty_mass[m, v] <- sum(PP[v, empty, m])
      }
      if (any(occupied)) {
        for (u in seq_len(K)) {
          B[v, u, m] <- sum(
            PP[v, occupied, m] *
              CountsV[u, occupied, m] /
              n_l[occupied]
          )
        }
      }

      if (isTRUE(conditional)) {
        ss <- sum(B[v, , m])
        if (ss > 0) B[v, , m] <- B[v, , m] / ss
      }
    }
  }

  B_mean <- apply(B, c(1L, 2L), mean)
  rownames(B_mean) <- colnames(B_mean) <- object$strata_levels
  colnames(empty_mass) <- object$strata_levels

  external_draws <- matrix(NA_real_, nrow = S, ncol = K)
  for (m in seq_len(S)) {
    external_draws[m, ] <- rowSums(B[, , m, drop = FALSE][, , 1L]) - diag(B[, , m])
  }
  colnames(external_draws) <- object$strata_levels

  out <- list(
    mean = B_mean,
    empty_mass_mean = colMeans(empty_mass),
    external_support_mean = colMeans(external_draws),
    conditional = isTRUE(conditional)
  )

  if (isTRUE(keep_draws)) {
    out$draws <- B
    out$empty_mass_draws <- empty_mass
    out$external_support_draws <- external_draws
  }

  class(out) <- c("strataLCM_support", "list")
  out
}

#' @export
print.strataLCM_support <- function(x, digits = 3L, ...) {
  cat("Posterior mean stratum support matrix\n")
  print(round(x$mean, digits = digits))
  cat("\nPosterior mean mass on globally empty components\n")
  print(round(x$empty_mass_mean, digits = digits))
  cat("\nPosterior mean external support (off-diagonal row sum)\n")
  print(round(x$external_support_mean, digits = digits))
  invisible(x)
}
