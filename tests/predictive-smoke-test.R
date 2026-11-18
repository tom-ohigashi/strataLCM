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
  n_draws = 30,
  seed = 2026
)

## One posterior draw -> data.frame
x1 <- lcm_draw_covariates(
  fit,
  stratum = "C",
  n = 100,
  draws = 3,
  seed = 9001,
  include_component = TRUE
)
stopifnot(is.data.frame(x1))
stopifnot(nrow(x1) == 100L)
stopifnot(identical(names(x1)[1:2], c("x1", "x2")))
stopifnot(all(x1$.component >= 1L & x1$.component <= 12L))
stopifnot(identical(attr(x1, "posterior_draw"), 3L))
stopifnot(identical(attr(x1, "stratum"), "C"))
stopifnot(all(is.finite(as.matrix(x1[, c("x1", "x2")]))) )

## Reproducibility
x1b <- lcm_draw_covariates(
  fit,
  stratum = "C",
  n = 100,
  draws = 3,
  seed = 9001,
  include_component = TRUE
)
stopifnot(identical(x1, x1b))

## Multiple posterior draws -> named list
xl <- lcm_draw_covariates(
  fit,
  stratum = "A",
  n = 20,
  draws = c(1, 4, 7),
  seed = 9002
)
stopifnot(is.list(xl), length(xl) == 3L)
stopifnot(identical(names(xl), c("draw_1", "draw_4", "draw_7")))
stopifnot(all(vapply(xl, nrow, integer(1)) == 20L))

## Fitted scale and original scale should be exactly related if the same RNG is used.
xs <- lcm_draw_covariates(
  fit,
  stratum = "B",
  n = 50,
  draws = 5,
  seed = 9003,
  original_scale = FALSE,
  include_component = TRUE
)
xo <- lcm_draw_covariates(
  fit,
  stratum = "B",
  n = 50,
  draws = 5,
  seed = 9003,
  original_scale = TRUE,
  include_component = TRUE
)
reconstructed <- sweep(as.matrix(xs[, c("x1", "x2")]), 2, fit$continuous$scale, "*")
reconstructed <- sweep(reconstructed, 2, fit$continuous$center, "+")
stopifnot(max(abs(reconstructed - as.matrix(xo[, c("x1", "x2")]))) < 1e-12)
stopifnot(identical(xs$.component, xo$.component))

## Mixed-data path
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
  n_draws = 20,
  seed = 2027
)

xm <- lcm_draw_covariates(
  fit_mixed,
  stratum = "C",
  n = 100,
  draws = 2,
  seed = 9004
)
stopifnot(is.factor(xm$sex), is.factor(xm$site))
stopifnot(identical(levels(xm$sex), fit_mixed$categorical$levels[[1]]))
stopifnot(identical(levels(xm$site), fit_mixed$categorical$levels[[2]]))
stopifnot(all(as.character(xm$sex) %in% c("0", "1")))
stopifnot(all(as.character(xm$site) %in% c("X", "Y", "Z")))

cat("strataLCM posterior-predictive smoke test passed.\n")
