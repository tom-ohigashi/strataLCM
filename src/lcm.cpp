// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

namespace {

inline bool is_integerish(double x, double tol = 1e-8) {
  return std::isfinite(x) && std::abs(x - std::round(x)) < tol;
}

void validate_common_inputs(const arma::mat& W_cont,
                            const arma::vec& v,
                            const arma::vec& initial_Z,
                            unsigned int mK,
                            double Alpha,
                            double mu_sd,
                            double prec_shape,
                            double prec_rate,
                            unsigned int NTHIN,
                            unsigned int NOUTSAMPLE) {
  if (W_cont.n_rows == 0 || W_cont.n_cols == 0) {
    Rcpp::stop("W_cont must have at least one row and one column.");
  }
  if (v.n_elem != W_cont.n_rows || initial_Z.n_elem != W_cont.n_rows) {
    Rcpp::stop("W_cont, v, and initial_Z must have the same number of observations.");
  }
  if (!W_cont.is_finite()) {
    Rcpp::stop("W_cont contains NA/NaN/Inf. Handle missing/non-finite values before fitting.");
  }
  if (mK < 1) Rcpp::stop("mK must be >= 1.");
  if (!(Alpha > 0.0) || !std::isfinite(Alpha)) Rcpp::stop("Alpha must be finite and > 0.");
  if (!(mu_sd > 0.0) || !std::isfinite(mu_sd)) Rcpp::stop("mu_sd must be finite and > 0.");
  if (!(prec_shape > 0.0) || !std::isfinite(prec_shape)) Rcpp::stop("prec_shape must be finite and > 0.");
  if (!(prec_rate > 0.0) || !std::isfinite(prec_rate)) Rcpp::stop("prec_rate must be finite and > 0.");
  if (NTHIN < 1) Rcpp::stop("NTHIN must be >= 1.");
  if (NOUTSAMPLE < 1) Rcpp::stop("NOUTSAMPLE must be >= 1.");

  // v must be 1,...,M without gaps because it is used as a row index.
  unsigned int vmax = 0;
  std::vector<bool> seen_v;
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    if (!is_integerish(v(i)) || v(i) < 1.0) {
      Rcpp::stop("v must contain positive integer labels 1,...,M.");
    }
    vmax = std::max(vmax, static_cast<unsigned int>(std::llround(v(i))));
  }
  seen_v.assign(vmax + 1, false);
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    seen_v[static_cast<unsigned int>(std::llround(v(i)))] = true;
  }
  for (unsigned int g = 1; g <= vmax; ++g) {
    if (!seen_v[g]) Rcpp::stop("v labels must be consecutive integers 1,...,M without gaps.");
  }

  for (arma::uword i = 0; i < initial_Z.n_elem; ++i) {
    if (!is_integerish(initial_Z(i)) || initial_Z(i) < 1.0 || initial_Z(i) > static_cast<double>(mK)) {
      Rcpp::stop("initial_Z must contain integer component labels in 1,...,mK.");
    }
  }
}

arma::uvec component_counts(const arma::uvec& z, unsigned int mK) {
  arma::uvec out(mK, arma::fill::zeros);
  for (arma::uword i = 0; i < z.n_elem; ++i) {
    if (z(i) >= mK) Rcpp::stop("Internal error: component label out of range.");
    out(z(i))++;
  }
  return out;
}

arma::vec rdirichlet_r(const arma::vec& alpha) {
  if (alpha.n_elem == 0) Rcpp::stop("Dirichlet parameter vector is empty.");
  for (arma::uword j = 0; j < alpha.n_elem; ++j) {
    if (!(alpha(j) > 0.0) || !std::isfinite(alpha(j))) {
      Rcpp::stop("All Dirichlet parameters must be finite and > 0.");
    }
  }

  arma::vec y(alpha.n_elem);
  // Very sparse Dirichlet parameters can make rgamma underflow to zero.
  // Retrying protects against the numerically degenerate all-zero draw.
  for (int attempt = 0; attempt < 100; ++attempt) {
    double total = 0.0;
    for (arma::uword j = 0; j < alpha.n_elem; ++j) {
      y(j) = R::rgamma(alpha(j), 1.0); // shape, scale
      total += y(j);
    }
    if (std::isfinite(total) && total > 0.0) return y / total;
  }
  Rcpp::stop("Dirichlet sampling failed numerically (all Gamma draws were zero/non-finite). Consider a less sparse Alpha.");
  return y; // unreachable
}

unsigned int sample_categorical_r(const arma::vec& prob) {
  if (prob.n_elem == 0) Rcpp::stop("Cannot sample from an empty probability vector.");
  double total = arma::accu(prob);
  if (!(total > 0.0) || !std::isfinite(total) || !prob.is_finite()) {
    Rcpp::stop("Invalid categorical probabilities encountered in cluster update.");
  }
  const double u = R::runif(0.0, 1.0);
  double cum = 0.0;
  for (arma::uword k = 0; k < prob.n_elem; ++k) {
    if (prob(k) < 0.0) Rcpp::stop("Negative categorical probability encountered.");
    cum += prob(k) / total;
    if (u <= cum || k + 1 == prob.n_elem) return static_cast<unsigned int>(k);
  }
  return static_cast<unsigned int>(prob.n_elem - 1); // unreachable safeguard
}

arma::uvec validate_and_get_cat_levels(const arma::mat& W_cat, arma::uword n_expected) {
  if (W_cat.n_rows != n_expected) {
    Rcpp::stop("W_cat and W_cont must have the same number of rows.");
  }
  if (W_cat.n_cols == 0) {
    Rcpp::stop("W_cat has zero columns. Use LCM_cpp1() for continuous-only data.");
  }
  if (!W_cat.is_finite()) {
    Rcpp::stop("W_cat contains NA/NaN/Inf. Handle missing values before fitting.");
  }

  arma::uvec levels(W_cat.n_cols);
  for (arma::uword j = 0; j < W_cat.n_cols; ++j) {
    unsigned int max_level = 0;
    for (arma::uword i = 0; i < W_cat.n_rows; ++i) {
      const double x = W_cat(i, j);
      if (!is_integerish(x) || x < 1.0) {
        Rcpp::stop("Each W_cat column must be coded as consecutive positive integers 1,...,C (not 0/1).");
      }
      max_level = std::max(max_level, static_cast<unsigned int>(std::llround(x)));
    }
    std::vector<bool> seen(max_level + 1, false);
    for (arma::uword i = 0; i < W_cat.n_rows; ++i) {
      seen[static_cast<unsigned int>(std::llround(W_cat(i, j)))] = true;
    }
    for (unsigned int c = 1; c <= max_level; ++c) {
      if (!seen[c]) {
        Rcpp::stop("Each W_cat column must use consecutive observed codes 1,...,C without gaps.");
      }
    }
    levels(j) = max_level;
  }
  return levels;
}

inline double safe_log_prob(double p) {
  if (p < 0.0 || !std::isfinite(p)) Rcpp::stop("Invalid component probability encountered.");
  if (p == 0.0) return -std::numeric_limits<double>::infinity();
  return std::log(p);
}

void print_progress(unsigned int iter, unsigned int NBURN, unsigned int NTHIN,
                    unsigned int NOUTSAMPLE, int print) {
  if (print != 1) return;
  const unsigned int total_post = NTHIN * NOUTSAMPLE;
  if (iter <= NBURN && NBURN > 0) {
    const unsigned int step = std::max(1u, NBURN / 10u);
    if (iter == NBURN || iter % step == 0) {
      Rcpp::Rcout << "Burn-in " << static_cast<int>(std::round(100.0 * iter / NBURN)) << "% completed\n";
    }
  } else if (iter > NBURN) {
    const unsigned int post_iter = iter - NBURN;
    const unsigned int step = std::max(1u, total_post / 10u);
    if (post_iter == total_post || post_iter % step == 0) {
      Rcpp::Rcout << "MCMC " << static_cast<int>(std::round(100.0 * post_iter / total_post)) << "% completed\n";
    }
  }
}

} // anonymous namespace

// Mixed continuous/categorical latent cluster mixture sampler.
// Categorical columns must be integer-coded 1,...,C_j.
// [[Rcpp::export]]
Rcpp::List LCM_cpp(const arma::mat& W_cont,
                   const arma::mat& W_cat,
                   const arma::vec& v,
                   const arma::vec& initial_Z,
                   const unsigned int mK = 7,
                   const double Alpha = 0.1,
                   const double mu_sd = 1000,
                   const double prec_shape = 0.1,
                   const double prec_rate = 0.1,
                   const unsigned int NBURN = 4000,
                   const unsigned int NTHIN = 10,
                   const unsigned int NOUTSAMPLE = 4000,
                   const int print = 0) {
  validate_common_inputs(W_cont, v, initial_Z, mK, Alpha, mu_sd, prec_shape, prec_rate, NTHIN, NOUTSAMPLE);
  const arma::uvec level_cat = validate_and_get_cat_levels(W_cat, W_cont.n_rows);

  const unsigned int p_cont = W_cont.n_cols;
  const unsigned int p_cat = W_cat.n_cols;
  const unsigned int M = static_cast<unsigned int>(arma::max(v));
  const unsigned int N = v.n_elem;
  const double alpha = Alpha / static_cast<double>(mK);
  const double base_pi = 1.0;
  const unsigned int level_max = level_cat.max();

  arma::uvec v0(N);
  arma::uvec Z(N);
  for (unsigned int i = 0; i < N; ++i) {
    v0(i) = static_cast<unsigned int>(std::llround(v(i))) - 1u;
    Z(i) = static_cast<unsigned int>(std::llround(initial_Z(i))) - 1u;
  }

  arma::mat Mu(p_cont, mK, arma::fill::zeros);
  arma::mat Sig(p_cont, mK, arma::fill::ones);
  arma::cube Pi(level_max, mK, p_cat, arma::fill::zeros);
  arma::mat PP(M, mK, arma::fill::ones);
  PP /= static_cast<double>(mK);

  // MCMC output
  arma::cube Mu_out(p_cont, mK, NOUTSAMPLE, arma::fill::zeros);
  arma::cube Sig_out(p_cont, mK, NOUTSAMPLE, arma::fill::zeros);
  std::vector<arma::cube> Pi_out;
  Pi_out.reserve(p_cat);
  for (unsigned int j = 0; j < p_cat; ++j) {
    Pi_out.emplace_back(level_max, mK, NOUTSAMPLE, arma::fill::zeros);
  }
  arma::cube PP_out(M, mK, NOUTSAMPLE, arma::fill::zeros);
  arma::umat Counts_out(mK, NOUTSAMPLE, arma::fill::zeros);
  // Stratum-by-component allocation counts for each saved posterior draw.
  // Dimensions in R: M strata x mK components x NOUTSAMPLE draws.
  arma::ucube CountsV_out(M, mK, NOUTSAMPLE, arma::fill::zeros);

  arma::vec log_dens(mK), prob(mK), dirpar, precision, post_prec, post_num;
  unsigned int count = 0;
  const unsigned int total_iter = NBURN + NTHIN * NOUTSAMPLE;

  for (unsigned int nit = 0; nit < total_iter; ++nit) {
    const arma::uvec counts = component_counts(Z, mK);

    // Update component-specific parameters.
    for (unsigned int k = 0; k < mK; ++k) {
      const arma::uword nk = counts(k);

      if (nk == 0) {
        // Empty component: draw from prior.
        for (unsigned int j = 0; j < p_cont; ++j) {
          Mu(j, k) = R::rnorm(0.0, mu_sd);
          const double lambda = R::rgamma(prec_shape, 1.0 / prec_rate); // precision
          if (!(lambda > 0.0) || !std::isfinite(lambda)) Rcpp::stop("Non-finite precision draw from prior.");
          Sig(j, k) = std::sqrt(1.0 / lambda);
        }
        for (unsigned int j = 0; j < p_cat; ++j) {
          dirpar = arma::ones<arma::vec>(level_cat(j)) * base_pi;
          Pi.slice(j).col(k).subvec(0, level_cat(j) - 1) = rdirichlet_r(dirpar);
        }
      } else {
        // Occupied component, INCLUDING singleton components (nk == 1).
        const arma::uvec idx = arma::find(Z == k);
        const arma::mat sW_cont = W_cont.rows(idx);
        const arma::mat sW_cat = W_cat.rows(idx);

        // Mu | sigma, data
        precision = 1.0 / arma::square(Sig.col(k));
        post_prec = static_cast<double>(nk) * precision + 1.0 / (mu_sd * mu_sd);
        post_num = arma::sum(sW_cont, 0).t() % precision;
        for (unsigned int j = 0; j < p_cont; ++j) {
          const double post_mean = post_num(j) / post_prec(j);
          const double post_sd = std::sqrt(1.0 / post_prec(j));
          Mu(j, k) = R::rnorm(post_mean, post_sd);
        }

        // sigma^2 | mu, data via precision ~ Gamma(shape, scale=1/rate)
        const arma::mat resid = sW_cont.each_row() - Mu.col(k).t();
        const arma::vec ss = arma::sum(arma::square(resid), 0).t();
        const double shape = prec_shape + 0.5 * static_cast<double>(nk);
        const arma::vec rate = prec_rate + 0.5 * ss;
        for (unsigned int j = 0; j < p_cont; ++j) {
          const double lambda = R::rgamma(shape, 1.0 / rate(j));
          if (!(lambda > 0.0) || !std::isfinite(lambda)) Rcpp::stop("Non-finite precision draw in component update.");
          Sig(j, k) = std::sqrt(1.0 / lambda);
        }

        // Categorical probabilities.
        for (unsigned int j = 0; j < p_cat; ++j) {
          dirpar = arma::ones<arma::vec>(level_cat(j)) * base_pi;
          for (arma::uword r = 0; r < sW_cat.n_rows; ++r) {
            const unsigned int lev = static_cast<unsigned int>(std::llround(sW_cat(r, j))) - 1u;
            dirpar(lev) += 1.0;
          }
          Pi.slice(j).col(k).subvec(0, level_cat(j) - 1) = rdirichlet_r(dirpar);
        }
      }
    }

    // Update latent component labels Z.
    for (unsigned int i = 0; i < N; ++i) {
      for (unsigned int k = 0; k < mK; ++k) {
        double lp = safe_log_prob(PP(v0(i), k));
        for (unsigned int j = 0; j < p_cont; ++j) {
          lp += R::dnorm(W_cont(i, j), Mu(j, k), Sig(j, k), true);
        }
        for (unsigned int j = 0; j < p_cat; ++j) {
          const unsigned int lev = static_cast<unsigned int>(std::llround(W_cat(i, j))) - 1u;
          // FIX: categorical contribution must be on the log scale.
          lp += safe_log_prob(Pi(lev, k, j));
        }
        log_dens(k) = lp;
      }
      const double max_lp = log_dens.max();
      if (!std::isfinite(max_lp)) Rcpp::stop("All component log-probabilities are non-finite for an observation.");
      prob = arma::exp(log_dens - max_lp);
      const double s = arma::accu(prob);
      if (!(s > 0.0) || !std::isfinite(s)) Rcpp::stop("Failed to normalize cluster probabilities.");
      prob /= s;
      Z(i) = sample_categorical_r(prob); // uses R RNG; reproducible under set.seed()
    }

    // Update stratum-specific mixture weights.
    for (unsigned int g = 0; g < M; ++g) {
      const arma::uvec zg = Z.elem(arma::find(v0 == g));
      arma::vec dpar = arma::conv_to<arma::vec>::from(component_counts(zg, mK));
      dpar += alpha;
      PP.row(g) = rdirichlet_r(dpar).t();
    }

    // Save thinned posterior draws.
    const unsigned int iter1 = nit + 1;
    if (iter1 > NBURN && ((iter1 - NBURN) % NTHIN == 0)) {
      if (count >= NOUTSAMPLE) Rcpp::stop("Internal output indexing error.");
      Mu_out.slice(count) = Mu;
      Sig_out.slice(count) = Sig;
      for (unsigned int j = 0; j < p_cat; ++j) Pi_out[j].slice(count) = Pi.slice(j);
      PP_out.slice(count) = PP;
      Counts_out.col(count) = component_counts(Z, mK);
      for (unsigned int g = 0; g < M; ++g) {
        const arma::uvec zg = Z.elem(arma::find(v0 == g));
        CountsV_out.slice(count).row(g) = component_counts(zg, mK).t();
      }
      ++count;
    }

    print_progress(iter1, NBURN, NTHIN, NOUTSAMPLE, print);
  }

  Rcpp::List Pi_out_list(p_cat);
  for (unsigned int j = 0; j < p_cat; ++j) Pi_out_list[j] = Pi_out[j];

  return Rcpp::List::create(
    Rcpp::_["Mu"] = Mu_out,
    Rcpp::_["Sig"] = Sig_out,
    Rcpp::_["Pi"] = Pi_out_list,
    Rcpp::_["PP"] = PP_out,
    Rcpp::_["Counts"] = Counts_out,
    Rcpp::_["CountsV"] = CountsV_out,
    Rcpp::_["Hyperparameters"] = Rcpp::List::create(
      Rcpp::_["Alpha"] = Alpha,
      Rcpp::_["alpha_component"] = alpha,
      Rcpp::_["mu_sd"] = mu_sd,
      Rcpp::_["prec_shape"] = prec_shape,
      Rcpp::_["prec_rate"] = prec_rate
    )
  );
}

// Continuous-only latent cluster mixture sampler.
// [[Rcpp::export]]
Rcpp::List LCM_cpp1(const arma::mat& W_cont,
                    const arma::vec& v,
                    const arma::vec& initial_Z,
                    const unsigned int mK = 7,
                    const double Alpha = 0.1,
                    const double mu_sd = 1000,
                    const double prec_shape = 0.1,
                    const double prec_rate = 0.1,
                    const unsigned int NBURN = 4000,
                    const unsigned int NTHIN = 10,
                    const unsigned int NOUTSAMPLE = 4000,
                    const int print = 0) {
  validate_common_inputs(W_cont, v, initial_Z, mK, Alpha, mu_sd, prec_shape, prec_rate, NTHIN, NOUTSAMPLE);

  const unsigned int p_cont = W_cont.n_cols;
  const unsigned int M = static_cast<unsigned int>(arma::max(v));
  const unsigned int N = v.n_elem;
  const double alpha = Alpha / static_cast<double>(mK);

  arma::uvec v0(N);
  arma::uvec Z(N);
  for (unsigned int i = 0; i < N; ++i) {
    v0(i) = static_cast<unsigned int>(std::llround(v(i))) - 1u;
    Z(i) = static_cast<unsigned int>(std::llround(initial_Z(i))) - 1u;
  }

  arma::mat Mu(p_cont, mK, arma::fill::zeros);
  arma::mat Sig(p_cont, mK, arma::fill::ones);
  arma::mat PP(M, mK, arma::fill::ones);
  PP /= static_cast<double>(mK);

  arma::cube Mu_out(p_cont, mK, NOUTSAMPLE, arma::fill::zeros);
  arma::cube Sig_out(p_cont, mK, NOUTSAMPLE, arma::fill::zeros);
  arma::cube PP_out(M, mK, NOUTSAMPLE, arma::fill::zeros);
  arma::umat Counts_out(mK, NOUTSAMPLE, arma::fill::zeros);
  // Stratum-by-component allocation counts for each saved posterior draw.
  // Dimensions in R: M strata x mK components x NOUTSAMPLE draws.
  arma::ucube CountsV_out(M, mK, NOUTSAMPLE, arma::fill::zeros);

  arma::vec log_dens(mK), prob(mK), precision, post_prec, post_num;
  unsigned int count = 0;
  const unsigned int total_iter = NBURN + NTHIN * NOUTSAMPLE;

  for (unsigned int nit = 0; nit < total_iter; ++nit) {
    const arma::uvec counts = component_counts(Z, mK);

    for (unsigned int k = 0; k < mK; ++k) {
      const arma::uword nk = counts(k);
      if (nk == 0) {
        for (unsigned int j = 0; j < p_cont; ++j) {
          Mu(j, k) = R::rnorm(0.0, mu_sd);
          const double lambda = R::rgamma(prec_shape, 1.0 / prec_rate);
          if (!(lambda > 0.0) || !std::isfinite(lambda)) Rcpp::stop("Non-finite precision draw from prior.");
          Sig(j, k) = std::sqrt(1.0 / lambda);
        }
      } else {
        // FIX: singleton components are updated from their one observation.
        const arma::uvec idx = arma::find(Z == k);
        const arma::mat sW_cont = W_cont.rows(idx);

        precision = 1.0 / arma::square(Sig.col(k));
        post_prec = static_cast<double>(nk) * precision + 1.0 / (mu_sd * mu_sd);
        post_num = arma::sum(sW_cont, 0).t() % precision;
        for (unsigned int j = 0; j < p_cont; ++j) {
          Mu(j, k) = R::rnorm(post_num(j) / post_prec(j), std::sqrt(1.0 / post_prec(j)));
        }

        const arma::mat resid = sW_cont.each_row() - Mu.col(k).t();
        const arma::vec ss = arma::sum(arma::square(resid), 0).t();
        const double shape = prec_shape + 0.5 * static_cast<double>(nk);
        const arma::vec rate = prec_rate + 0.5 * ss;
        for (unsigned int j = 0; j < p_cont; ++j) {
          const double lambda = R::rgamma(shape, 1.0 / rate(j));
          if (!(lambda > 0.0) || !std::isfinite(lambda)) Rcpp::stop("Non-finite precision draw in component update.");
          Sig(j, k) = std::sqrt(1.0 / lambda);
        }
      }
    }

    for (unsigned int i = 0; i < N; ++i) {
      for (unsigned int k = 0; k < mK; ++k) {
        double lp = safe_log_prob(PP(v0(i), k));
        for (unsigned int j = 0; j < p_cont; ++j) {
          lp += R::dnorm(W_cont(i, j), Mu(j, k), Sig(j, k), true);
        }
        log_dens(k) = lp;
      }
      const double max_lp = log_dens.max();
      if (!std::isfinite(max_lp)) Rcpp::stop("All component log-probabilities are non-finite for an observation.");
      prob = arma::exp(log_dens - max_lp);
      const double s = arma::accu(prob);
      if (!(s > 0.0) || !std::isfinite(s)) Rcpp::stop("Failed to normalize cluster probabilities.");
      prob /= s;
      Z(i) = sample_categorical_r(prob);
    }

    for (unsigned int g = 0; g < M; ++g) {
      const arma::uvec zg = Z.elem(arma::find(v0 == g));
      arma::vec dpar = arma::conv_to<arma::vec>::from(component_counts(zg, mK));
      dpar += alpha;
      PP.row(g) = rdirichlet_r(dpar).t();
    }

    const unsigned int iter1 = nit + 1;
    if (iter1 > NBURN && ((iter1 - NBURN) % NTHIN == 0)) {
      if (count >= NOUTSAMPLE) Rcpp::stop("Internal output indexing error.");
      Mu_out.slice(count) = Mu;
      Sig_out.slice(count) = Sig;
      PP_out.slice(count) = PP;
      Counts_out.col(count) = component_counts(Z, mK);
      for (unsigned int g = 0; g < M; ++g) {
        const arma::uvec zg = Z.elem(arma::find(v0 == g));
        CountsV_out.slice(count).row(g) = component_counts(zg, mK).t();
      }
      ++count;
    }

    print_progress(iter1, NBURN, NTHIN, NOUTSAMPLE, print);
  }

  return Rcpp::List::create(
    Rcpp::_["Mu"] = Mu_out,
    Rcpp::_["Sig"] = Sig_out,
    Rcpp::_["PP"] = PP_out,
    Rcpp::_["Counts"] = Counts_out,
    Rcpp::_["CountsV"] = CountsV_out,
    Rcpp::_["Hyperparameters"] = Rcpp::List::create(
      Rcpp::_["Alpha"] = Alpha,
      Rcpp::_["alpha_component"] = alpha,
      Rcpp::_["mu_sd"] = mu_sd,
      Rcpp::_["prec_shape"] = prec_shape,
      Rcpp::_["prec_rate"] = prec_rate
    )
  );
}
