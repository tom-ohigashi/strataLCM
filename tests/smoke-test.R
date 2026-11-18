# Minimal smoke tests for the C++ core.
library(strataLCM)

set.seed(20260920)
n <- 40
v <- rep(1:2, each = n / 2)
W <- cbind(rnorm(n, rep(c(-1, 1), each = n / 2), 1), rnorm(n))
initial_Z <- rep(1:4, length.out = n)

set.seed(123)
f1 <- LCM_cpp1(W, v, initial_Z, mK = 4, Alpha = 0.1,
               mu_sd = 5, prec_shape = 1, prec_rate = 1, NBURN = 5, NTHIN = 1,
               NOUTSAMPLE = 3, print = 0)
set.seed(123)
f2 <- LCM_cpp1(W, v, initial_Z, mK = 4, Alpha = 0.1,
               mu_sd = 5, prec_shape = 1, prec_rate = 1, NBURN = 5, NTHIN = 1,
               NOUTSAMPLE = 3, print = 0)

stopifnot(identical(f1$Mu, f2$Mu))
stopifnot(identical(f1$Sig, f2$Sig))
stopifnot(identical(f1$PP, f2$PP))
stopifnot(all(colSums(f1$Counts) == n))

# Mixed-data smoke test. Note binary variable is coded 1/2, not 0/1.
Wcat <- cbind(ifelse(W[, 1] > 0, 2, 1), rep(1:2, length.out = n))
set.seed(456)
fm <- LCM_cpp(W, Wcat, v, initial_Z, mK = 4, Alpha = 0.1,
              mu_sd = 5, prec_shape = 1, prec_rate = 1, NBURN = 5, NTHIN = 1,
              NOUTSAMPLE = 3, print = 0)
stopifnot(all(colSums(fm$Counts) == n))
stopifnot(length(fm$Pi) == ncol(Wcat))
