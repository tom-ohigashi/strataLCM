library(strataLCM)

set.seed(101)
nv <- c(30, 12, 5)
strata <- rep(c("A", "B", "C"), nv)
W <- rbind(
  cbind(rnorm(nv[1], 0, 1), rnorm(nv[1], 0, 1)),
  cbind(rnorm(nv[2], 2, 1), rnorm(nv[2], 2, 1)),
  cbind(rnorm(nv[3], -2, 1), rnorm(nv[3], -2, 1))
)
colnames(W) <- c("x1", "x2")

fit <- lcm_fit(
  W_cont = W,
  strata = strata,
  n_components = 12,
  nburn = 50,
  thin = 1,
  n_draws = 50,
  seed = 2026
)

stopifnot(inherits(fit, "strataLCM_fit"))
stopifnot(identical(dim(fit$PP), c(3L, 12L, 50L)))
stopifnot(identical(dim(fit$CountsV), c(3L, 12L, 50L)))
stopifnot(all(apply(fit$CountsV, c(2, 3), sum) == fit$Counts))
stopifnot(all(apply(fit$CountsV, c(1, 3), sum)[1, ] == nv[1]))
stopifnot(all(apply(fit$CountsV, c(1, 3), sum)[2, ] == nv[2]))
stopifnot(all(apply(fit$CountsV, c(1, 3), sum)[3, ] == nv[3]))
stopifnot(max(abs(colMeans(scale(W)))) < 1e-10)

sup <- lcm_support(fit)
stopifnot(identical(dim(sup$mean), c(3L, 3L)))
stopifnot(all(is.finite(sup$mean)))
stopifnot(all(sup$external_support_mean >= 0))

# Mixed-data path: 0/1 and character categories are encoded internally.
Wcat <- data.frame(
  sex = rep(c(0, 1), length.out = nrow(W)),
  site = rep(c("X", "Y", "Z"), length.out = nrow(W))
)

fit_mixed <- lcm_fit(
  W_cont = W,
  W_cat = Wcat,
  strata = strata,
  n_components = 12,
  nburn = 25,
  thin = 1,
  n_draws = 25,
  seed = 2027
)

stopifnot(inherits(fit_mixed, "strataLCM_fit"))
stopifnot(length(fit_mixed$categorical$levels) == 2L)
stopifnot(identical(dim(fit_mixed$CountsV), c(3L, 12L, 25L)))

cat("strataLCM R-wrapper smoke test passed.\n")
