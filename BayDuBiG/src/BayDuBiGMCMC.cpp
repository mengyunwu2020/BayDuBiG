// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// OpenMP is optional.  When the file is compiled without -fopenmp (e.g. the
// default Apple clang toolchain), _OPENMP is undefined and the code falls back
// to a fully serial implementation that is numerically identical.
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// ============================================================================
// 1. Stochastic trace estimation of powers of the spatial weight matrix
//    (Hutchinson estimator). computeSingleTraceEstimate(W, i, v) estimates
//    tr(W^i) using a single random vector v.
// ============================================================================
double computeSingleTraceEstimate(const arma::sp_mat &W, int i, const arma::vec &v) {
  arma::vec Wv = v;
  for (int k = 0; k < i; ++k) {
    Wv = W * Wv;
  }
  return arma::dot(v, Wv);
}

// Approximate log|I - rho*W| = -sum_i rho^i * tr(W^i) / i using the precomputed traces.
double approximateLogDetPrecomputed(const std::vector<std::vector<double>>& precomputedTraces,
                                    double rho, int order, int m) {
  double logDetEstimate = 0.0;
  for (int i = 1; i <= order; ++i) {
    double traceEstimateSumForOrder = 0.0;
    for (int j = 0; j < m; ++j) {
      traceEstimateSumForOrder += precomputedTraces[i - 1][j];
    }
    logDetEstimate -= std::pow(rho, i) * traceEstimateSumForOrder / i;
  }
  return logDetEstimate / m;
}

// (Approximate) log-likelihood for each gene.
arma::vec log_likelihood(const arma::mat& y_norm, const arma::vec &rho, const arma::vec &sigma,
                         const arma::mat &Z, const arma::mat &psi, const arma::mat &WY,
                         const int &p, const int &n, const int &m, const int &order,
                         const std::vector<std::vector<double>>& precomputedTraces) {
  arma::vec log_likes(p);
  for (arma::uword j = 0; j < y_norm.n_cols; ++j) {
    double det_part = approximateLogDetPrecomputed(precomputedTraces, rho(j), order, m);
    double sgm_part = -(n / 2.0) * std::log(2.0 * M_PI * sigma(j));

    arma::vec e_part = y_norm.col(j) - rho(j) * WY.col(j) - Z * psi.col(j);
    double exp_part = -arma::dot(e_part, e_part) / (2.0 * sigma(j));
    log_likes(j) = det_part + sgm_part + exp_part;
  }
  return log_likes;
}

// ============================================================================
// 2. Group effect under the union semantics
//
//    A gene's group effect is zeta = 1 - prod_k (1 - tau[group_k]), so it is
//    active when at least one of the gene's member groups is selected.
//
//    union_tau      : computes a gene's zeta.
//    union_tau_excl : computes the group effect of the remaining groups after
//                     excluding the group g currently being updated; used in the
//                     Gibbs update of a single tau_g to compare the likelihood
//                     under tau_g = 1 versus tau_g = 0.
// ============================================================================
static inline double union_tau(const Rcpp::IntegerVector& groups, const arma::ivec& tau) {
  if (groups.length() == 0) return 0.0;
  double p = 1.0;
  for (int k = 0; k < groups.length(); ++k) {
    p *= (1.0 - static_cast<double>(tau[groups[k]]));
  }
  return 1.0 - p;
}

static inline double union_tau_excl(const std::vector<int>& groups, const arma::ivec& tau,
                                    int excl, int G) {
  double p = 1.0;
  bool any = false;
  for (const int& g : groups) {
    if (g != excl && g >= 0 && g < G) {
      any = true;
      p *= (1.0 - static_cast<double>(tau[g]));
    }
  }
  return any ? (1.0 - p) : 0.0;
}

// ============================================================================
// 3. Core MCMC sampler (under a single covariate basis Z)
// ============================================================================
Rcpp::List run_mcmc(const arma::mat& Y,
                    int n,
                    int p,
                    const Rcpp::List& gene_group_list,
                    int G,
                    const arma::mat& Z,
                    const arma::mat& WY,
                    int iter,
                    int burn,
                    double a_rho,
                    double b_rho,
                    double a_tau,
                    double b_tau,
                    double a_gamma,
                    double b_gamma,
                    double a_sigma,
                    double b_sigma,
                    const std::vector<std::vector<double>>& precomputedTraces,
                    int order,
                    int m,
                    double time_build_WY_ms,
                    const std::string& model_name,
                    bool informative_gamma = false,
                    double a_gamma_on = 10.0,
                    double b_gamma_on = 1.0,
                    double a_gamma_off = 1.0,
                    double b_gamma_off = 10.0) {

  int j, it;
  int z_dim = Z.n_cols;
  double hastings;
  arma::vec rho(p);
  arma::Col<int> tau(G);
  arma::Col<int> gamma(p);

  arma::mat psi(z_dim, p);
  arma::vec sigma(p);

  arma::vec rho_sum(p, arma::fill::zeros);
  arma::vec sigma_sum(p, arma::fill::zeros);
  std::vector<arma::vec> psi_sum(p, arma::vec(z_dim).fill(0));

  int effective_iterations = iter - burn;

  arma::vec loglike_sum(p, arma::fill::zeros);
  arma::vec loglike_gene_mean(p);

  arma::ivec tau_1(G, fill::zeros);
  arma::ivec tau_2(G, fill::zeros);
  arma::vec tau1_product(p, fill::zeros);   // zeta_{1j}
  arma::vec tau2_product(p, fill::zeros);   // zeta_{2j}
  arma::ivec gamma_1(p, fill::zeros);
  arma::ivec gamma_2(p, fill::zeros);

  arma::vec tau_gamma_1_sum(p, fill::zeros);
  arma::vec tau_gamma_2_sum(p, fill::zeros);

  double sigma_start = 0.01;
  std::vector<double> c(p, 0.1);
  std::vector<int> accept_count(p, 0);
  std::vector<int> total_count(p, 0);

  for (j = 0; j < p; j++) {
    rho(j) = 0;
    gamma_1(j) = 1;
    gamma_2(j) = 1;
    sigma(j) = sigma_start;
    for (int z = 0; z < z_dim; z++) {
      psi(z, j) = 1;
    }
  }

  for (int g = 0; g < G; g++) {
    tau_1(g) = 1;
    tau_2(g) = 1;
  }

  // Initialize zeta = union_tau(groups)
  for (int j = 0; j < p; j++) {
    Rcpp::IntegerVector gene_groups = gene_group_list[j];
    tau1_product(j) = union_tau(gene_groups, tau_1);
    tau2_product(j) = union_tau(gene_groups, tau_2);
  }

  arma::mat y_norm = Y;

  // Cholesky factorization of Z'Z (ridge-regularized to guarantee positive
  // definiteness); reused in the psi updates below.
  arma::mat Zt = Z.t();
  arma::mat ZtZ = Zt * Z;
  const double ridge_eps = 1e-8;
  ZtZ.diag() += ridge_eps;
  arma::mat Rchol;
  bool chol_ok = arma::chol(Rchol, ZtZ, "upper");
  if (!chol_ok) {
    ZtZ.diag() += 1e-6;
    chol_ok = arma::chol(Rchol, ZtZ, "upper");
    if (!chol_ok) {
      stop("Cholesky factorization failed for ZtZ even after ridge regularization.");
    }
  }

  // The first base_cols columns are the baseline covariates (XWX); the last two
  // columns are the spatial-coordinate terms. When XWX is absent, Z has only the
  // two coordinate columns and base_cols = 0.
  int base_cols = (z_dim > 2) ? (z_dim - 2) : 0;
  arma::mat Z_base;
  if (base_cols > 0) {
    Z_base = Z.cols(0, base_cols - 1);
  } else {
    Z_base = arma::mat(Z.n_rows, 0);
  }

  // Precompute the groups each gene belongs to, and the genes within each group.
  std::vector<std::vector<int>> group_genes(G);
  std::vector<std::vector<int>> gene_groups_vec(p);
  for (int j = 0; j < p; ++j) {
    if (gene_group_list.length() > j && !Rf_isNull(gene_group_list[j])) {
      Rcpp::IntegerVector groups = gene_group_list[j];
      gene_groups_vec[j] = Rcpp::as<std::vector<int>>(groups);
      for (int group : groups) {
        if (group >= 0 && group < G) {
          group_genes[group].push_back(j);
        }
      }
    }
  }

  Rcout << "[Model] " << model_name;
  int checkpoint_iter = std::max(1, iter / 5);

  // ----------------------------------------------------------------------
  // Main MCMC loop
  // ----------------------------------------------------------------------
  for (it = 0; it < iter; it++) {
    // Precompute tildeY = Y - rho * WY and the baseline part R_base = Z_base * psi_base.
    arma::mat tildeY_all(y_norm.n_rows, p);
    arma::mat R_base_all(y_norm.n_rows, p);
    for (int j_pre = 0; j_pre < p; ++j_pre) {
      tildeY_all.col(j_pre) = y_norm.col(j_pre) - rho(j_pre) * WY.col(j_pre);
      if (base_cols > 0) {
        R_base_all.col(j_pre) = Z_base * psi.col(j_pre).rows(0, base_cols - 1);
      } else {
        R_base_all.col(j_pre).zeros();
      }
    }

    // Recompute each gene's zeta = union_tau(groups).
    for (int j = 0; j < p; ++j) {
      Rcpp::IntegerVector gene_groups = gene_group_list[j];
      tau1_product(j) = union_tau(gene_groups, tau_1);
      tau2_product(j) = union_tau(gene_groups, tau_2);
    }

    // --- Update tau_1 (group-level indicator for the first spatial trend) ---
    std::vector<double> posterior_probs_tau1(G);
    for (int g_idx = 0; g_idx < G; ++g_idx) {
      double pi_k1 = double(a_tau) / (a_tau + b_tau);
      double log_prob_tau_k1_1 = std::log(pi_k1);
      double log_prob_tau_k1_0 = std::log(1.0 - pi_k1);
      double log_likelihood_tau_k1_1 = 0.0;
      double log_likelihood_tau_k1_0 = 0.0;

      const std::vector<int>& genes_in_group = group_genes[g_idx];
#ifdef _OPENMP
#pragma omp parallel for reduction(+:log_likelihood_tau_k1_1, log_likelihood_tau_k1_0)
#endif
      for (size_t idx = 0; idx < genes_in_group.size(); ++idx) {
        int j = genes_in_group[idx];
        const std::vector<int>& gene_groups = gene_groups_vec[j];
        double other_groups_effect = union_tau_excl(gene_groups, tau_1, g_idx, G);

        arma::vec tilde_y_j = tildeY_all.col(j);

        // When tau_1[g_idx] = 1, zeta is always 1.
        arma::vec e_j_1 = tilde_y_j - R_base_all.col(j)
          - gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
          - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);

        // When tau_1[g_idx] = 0, zeta = the union effect of the remaining groups.
        arma::vec e_j_0 = tilde_y_j - R_base_all.col(j)
          - other_groups_effect * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
          - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);

        log_likelihood_tau_k1_1 += -arma::dot(e_j_1, e_j_1) / (2.0 * sigma(j));
        log_likelihood_tau_k1_0 += -arma::dot(e_j_0, e_j_0) / (2.0 * sigma(j));
      }

      log_prob_tau_k1_1 += log_likelihood_tau_k1_1;
      log_prob_tau_k1_0 += log_likelihood_tau_k1_0;

      double max_log_prob = std::max(log_prob_tau_k1_1, log_prob_tau_k1_0);
      double prob_tau_k1_1 = std::exp(log_prob_tau_k1_1 - max_log_prob);
      double prob_tau_k1_0 = std::exp(log_prob_tau_k1_0 - max_log_prob);

      posterior_probs_tau1[g_idx] = prob_tau_k1_1 / (prob_tau_k1_0 + prob_tau_k1_1);
      tau_1(g_idx) = R::rbinom(1, posterior_probs_tau1[g_idx]);
    }

    // --- Update tau_2 (group-level indicator for the second spatial trend) ---
    std::vector<double> posterior_probs_tau2(G);
    for (int g_idx = 0; g_idx < G; ++g_idx) {
      double pi_k2 = double(a_tau) / (a_tau + b_tau);
      double log_prob_tau_k2_1 = std::log(pi_k2);
      double log_prob_tau_k2_0 = std::log(1.0 - pi_k2);
      double log_likelihood_tau_k2_1 = 0.0;
      double log_likelihood_tau_k2_0 = 0.0;

      const std::vector<int>& genes_in_group = group_genes[g_idx];
#ifdef _OPENMP
#pragma omp parallel for reduction(+:log_likelihood_tau_k2_1, log_likelihood_tau_k2_0)
#endif
      for (size_t idx = 0; idx < genes_in_group.size(); ++idx) {
        int j = genes_in_group[idx];
        const std::vector<int>& gene_groups = gene_groups_vec[j];
        double other_groups_effect = union_tau_excl(gene_groups, tau_2, g_idx, G);

        arma::vec tilde_y_j = tildeY_all.col(j);

        arma::vec e_j_1 = tilde_y_j - R_base_all.col(j)
          - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
          - gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);

        arma::vec e_j_0 = tilde_y_j - R_base_all.col(j)
          - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
          - other_groups_effect * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);

        log_likelihood_tau_k2_1 += -arma::dot(e_j_1, e_j_1) / (2.0 * sigma(j));
        log_likelihood_tau_k2_0 += -arma::dot(e_j_0, e_j_0) / (2.0 * sigma(j));
      }

      log_prob_tau_k2_1 += log_likelihood_tau_k2_1;
      log_prob_tau_k2_0 += log_likelihood_tau_k2_0;

      double max_log_prob = std::max(log_prob_tau_k2_1, log_prob_tau_k2_0);
      double prob_tau_k2_1 = std::exp(log_prob_tau_k2_1 - max_log_prob);
      double prob_tau_k2_0 = std::exp(log_prob_tau_k2_0 - max_log_prob);

      posterior_probs_tau2[g_idx] = prob_tau_k2_1 / (prob_tau_k2_0 + prob_tau_k2_1);
      tau_2(g_idx) = R::rbinom(1, posterior_probs_tau2[g_idx]);
    }

    // --- Update gamma_1 (gene level) ---
    for (int j = 0; j < p; ++j) {
      double a_gamma_use = a_gamma;
      double b_gamma_use = b_gamma;
      if (informative_gamma) {
        if (tau1_product(j) > 0.5) {
          a_gamma_use = a_gamma_on;
          b_gamma_use = b_gamma_on;
        } else {
          a_gamma_use = a_gamma_off;
          b_gamma_use = b_gamma_off;
        }
      }
      double pi_j1 = double(a_gamma_use) / (a_gamma_use + b_gamma_use);
      double log_prob_gamma_j1_1 = std::log(pi_j1);
      double log_prob_gamma_j1_0 = std::log(1.0 - pi_j1);
      arma::vec tilde_y_jk = tildeY_all.col(j);
      arma::vec e_j1_1 = tilde_y_jk - R_base_all.col(j)
        - tau1_product(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
        - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
      arma::vec e_j1_0 = tilde_y_jk - R_base_all.col(j)
        - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
      log_prob_gamma_j1_1 += -arma::dot(e_j1_1, e_j1_1) / (2.0 * sigma(j));
      log_prob_gamma_j1_0 += -arma::dot(e_j1_0, e_j1_0) / (2.0 * sigma(j));
      double max_log_prob = std::max(log_prob_gamma_j1_1, log_prob_gamma_j1_0);
      double prob_gamma_j1_1 = std::exp(log_prob_gamma_j1_1 - max_log_prob);
      double prob_gamma_j1_0 = std::exp(log_prob_gamma_j1_0 - max_log_prob);
      double posterior_prob_gamma_j1_1 = prob_gamma_j1_1 / (prob_gamma_j1_1 + prob_gamma_j1_0);
      gamma_1(j) = R::rbinom(1, posterior_prob_gamma_j1_1);
    }

    // --- Update gamma_2 (gene level) ---
    for (int j = 0; j < p; ++j) {
      double a_gamma_use = a_gamma;
      double b_gamma_use = b_gamma;
      if (informative_gamma) {
        if (tau2_product(j) > 0.5) {
          a_gamma_use = a_gamma_on;
          b_gamma_use = b_gamma_on;
        } else {
          a_gamma_use = a_gamma_off;
          b_gamma_use = b_gamma_off;
        }
      }
      double pi_j2 = double(a_gamma_use) / (a_gamma_use + b_gamma_use);
      double log_prob_gamma_j2_1 = std::log(pi_j2);
      double log_prob_gamma_j2_0 = std::log(1.0 - pi_j2);
      arma::vec tilde_y_jk = tildeY_all.col(j);
      arma::vec e_j2_1 = tilde_y_jk - R_base_all.col(j)
        - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
        - tau2_product(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
      arma::vec e_j2_0 = tilde_y_jk - R_base_all.col(j)
        - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j);
      log_prob_gamma_j2_1 += -arma::dot(e_j2_1, e_j2_1) / (2.0 * sigma(j));
      log_prob_gamma_j2_0 += -arma::dot(e_j2_0, e_j2_0) / (2.0 * sigma(j));
      double max_log_prob = std::max(log_prob_gamma_j2_1, log_prob_gamma_j2_0);
      double prob_gamma_j2_1 = std::exp(log_prob_gamma_j2_1 - max_log_prob);
      double prob_gamma_j2_0 = std::exp(log_prob_gamma_j2_0 - max_log_prob);
      double posterior_prob_gamma_j2_1 = prob_gamma_j2_1 / (prob_gamma_j2_1 + prob_gamma_j2_0);
      gamma_2(j) = R::rbinom(1, posterior_prob_gamma_j2_1);
    }

    // --- Update psi (covariate coefficients via Cholesky; serial to keep RNG correct) ---
    for (int j = 0; j < p; j++) {
      arma::vec y_norm_tilde = y_norm.col(j) - rho(j) * WY.col(j);
      arma::vec b = Zt * y_norm_tilde;
      arma::vec u = arma::solve(arma::trimatl(Rchol.t()), b, arma::solve_opts::fast);
      arma::vec r = arma::solve(arma::trimatu(Rchol), u, arma::solve_opts::fast);

      r(z_dim - 2) *= tau1_product(j) * gamma_1(j);
      r(z_dim - 1) *= tau2_product(j) * gamma_2(j);

      arma::vec z = arma::randn<arma::vec>(z_dim);
      arma::vec a = arma::solve(arma::trimatl(Rchol.t()), z, arma::solve_opts::fast);
      arma::vec w = arma::solve(arma::trimatu(Rchol), a, arma::solve_opts::fast);
      psi.col(j) = r + std::sqrt(sigma(j)) * w;

      if (tau1_product(j) * gamma_1(j) < 1e-10) {
        psi(z_dim - 2, j) = 0.0;
      }
      if (tau2_product(j) * gamma_2(j) < 1e-10) {
        psi(z_dim - 1, j) = 0.0;
      }
    }

    // --- Update sigma (inverse Gamma) ---
    for (int j = 0; j < p; j++) {
      arma::vec psi_corrected = psi.col(j);
      psi_corrected(z_dim - 2) *= tau1_product(j) * gamma_1(j);
      psi_corrected(z_dim - 1) *= tau2_product(j) * gamma_2(j);
      arma::vec Z_psi_part = Z * psi_corrected;

      arma::vec e_kj = y_norm.col(j) - rho(j) * WY.col(j) - Z_psi_part;
      if (e_kj.has_nan() || e_kj.has_inf()) {
        Rcout << "e_kj has invalid values for j: " << j << std::endl;
        continue;
      }
      double alpha = n / 2.0 + a_sigma / 2.0;
      double beta = arma::dot(e_kj, e_kj) / 2.0 + b_sigma / 2.0;
      sigma(j) = 1.0 / R::rgamma(alpha, 1.0 / beta);
    }

    // --- Update rho (random-walk Metropolis) ---
    for (int j = 0; j < p; j++) {
      hastings = 0;
      double rho_temp = rho(j) + c[j] * R::rnorm(0, 1);
      double det_current = approximateLogDetPrecomputed(precomputedTraces, rho(j), order, m);
      double det_proposed = approximateLogDetPrecomputed(precomputedTraces, rho_temp, order, m);
      double det_ratio = std::exp(det_proposed - det_current);
      arma::vec e_kj_old = y_norm.col(j) - Z * psi.col(j) - rho(j) * WY.col(j);
      arma::vec e_kj_new = y_norm.col(j) - Z * psi.col(j) - rho_temp * WY.col(j);
      double error_ratio = std::exp(-arma::dot(e_kj_new, e_kj_new) / (2.0 * sigma(j))
                                     + arma::dot(e_kj_old, e_kj_old) / (2.0 * sigma(j)));
      hastings += error_ratio * det_ratio;
      if (hastings >= double(rand() % 10001) / 10000.0 && rho_temp > -1 && rho_temp < 1) {
        rho(j) = rho_temp;
        accept_count[j]++;
      }
      total_count[j]++;
      double accept_rate = double(accept_count[j]) / total_count[j];
      if (accept_rate > 0.30) {
        c[j] *= 1.1;
      } else if (accept_rate < 0.20) {
        c[j] /= 1.1;
      }
    }

    if ((it + 1) % checkpoint_iter == 0 || (it + 1) == iter) {
      int pct = (it + 1) * 100 / iter;
      Rcout << "[Model] " << model_name << " progress: " << pct << "%\n";
    }

    if (it >= burn) {
      for (j = 0; j < p; j++) {
        rho_sum(j) += rho(j);
        sigma_sum(j) += sigma(j);
        psi_sum[j] += psi.col(j);

        tau_gamma_1_sum(j) += tau1_product(j) * gamma_1(j);
        tau_gamma_2_sum(j) += tau2_product(j) * gamma_2(j);
      }
      arma::vec gene_loglike = log_likelihood(y_norm, rho, sigma, Z, psi, WY,
                                              p, n, m, order, precomputedTraces);
      for (int j = 0; j < p; j++) {
        loglike_sum(j) += gene_loglike(j);
      }
    }
  }

  arma::vec rho_mean = rho_sum / static_cast<double>(effective_iterations);
  arma::vec sigma_mean = sigma_sum / static_cast<double>(effective_iterations);
  for (int j = 0; j < p; j++) {
    psi_sum[j] /= static_cast<double>(effective_iterations);
  }

  arma::vec tau_gamma_1_mean = tau_gamma_1_sum / static_cast<double>(effective_iterations);
  arma::vec tau_gamma_2_mean = tau_gamma_2_sum / static_cast<double>(effective_iterations);

  // PPI: take the larger of the two group-effect posterior probabilities
  // (equivalent to the mean of max(zeta1*gamma1, zeta2*gamma2)).
  arma::vec tau_gamma(p, arma::fill::zeros);
  for (int i = 0; i < p; ++i) {
    tau_gamma(i) = std::max(tau_gamma_1_mean(i), tau_gamma_2_mean(i));
  }

  double total_loglike_sum = 0.0;
  for (int j = 0; j < p; j++) {
    loglike_gene_mean(j) = loglike_sum(j) / static_cast<double>(effective_iterations);
    total_loglike_sum += loglike_sum(j);
  }
  double total_loglike_mean = total_loglike_sum / static_cast<double>(effective_iterations);

  return Rcpp::List::create(
    Rcpp::Named("tau_gamma_upper") = tau_gamma,
    Rcpp::Named("gene_log_likelihood_mean") = loglike_gene_mean,
    Rcpp::Named("avg_log_likelihood") = total_loglike_mean,
    Rcpp::Named("rho_mean") = rho_mean,
    Rcpp::Named("sigma_mean") = sigma_mean,
    Rcpp::Named("psi_mean") = psi_sum
  );
}

// Safe standardization: return a zero vector when the standard deviation is near
// zero, to avoid NaN from 0/0.
static inline arma::vec safe_standardize(const arma::vec& x) {
  double mu = arma::mean(x);
  double sd = arma::stddev(x);
  if (!std::isfinite(sd) || sd < 1e-12) {
    return arma::vec(x.n_elem, arma::fill::zeros);
  }
  return (x - mu) / sd;
}

// ============================================================================
// 4. Exported entry point: BayDuBiG_mcmc
// ============================================================================
// [[Rcpp::export]]
Rcpp::List BayDuBiG_mcmc(const arma::mat& Y,
                         const Rcpp::List& gene_group_list,
                         const Rcpp::Nullable<Rcpp::NumericMatrix>& XWX,
                         const arma::vec& X_coord,
                         const arma::vec& Y_coord,
                         const arma::sp_mat& weights,
                         int iter,
                         int burn,
                         double a_rho = 1.00,
                         double b_rho = 1.00,
                         double a_tau = 1,
                         double b_tau = 1,
                         double a_gamma = 1,
                         double b_gamma = 1,
                         double a_sigma = 0.01,
                         double b_sigma = 0.01,
                         bool informative_gamma = false,
                         double a_gamma_on = 10.0,
                         double b_gamma_on = 1.0,
                         double a_gamma_off = 1.0,
                         double b_gamma_off = 10.0) {

#ifdef _OPENMP
  Rcout << "OpenMP is enabled. Number of threads: " << omp_get_max_threads() << std::endl;
#endif

  int n = Y.n_rows;
  int p = Y.n_cols;

  int order = 5;
  int m = 1;

  // Precompute stochastic trace estimates of tr(W^i).
  std::vector<std::vector<double>> precomputedTraces(order, std::vector<double>(m));
  for (int j = 0; j < m; ++j) {
    arma::vec v = arma::randn<arma::vec>(n);
    for (int i = 0; i < order; ++i) {
      precomputedTraces[i][j] = computeSingleTraceEstimate(weights, i + 1, v);
    }
  }

  arma::mat Z_original;
  arma::mat Z_periodic;
  arma::mat Z_exponential;

  if (XWX.isNotNull()) {
    Rcpp::NumericMatrix XWX_R = XWX.get();
    arma::mat XWX_mat(XWX_R.begin(), XWX_R.nrow(), XWX_R.ncol(), false);

    Z_original = arma::join_horiz(XWX_mat, arma::join_horiz(X_coord, Y_coord));
    Z_periodic = arma::join_horiz(XWX_mat,
                                  arma::join_horiz(cos(2 * M_PI * X_coord),
                                                   cos(2 * M_PI * Y_coord)));

    arma::vec exp_x2 = safe_standardize(exp(-X_coord % X_coord));
    arma::vec exp_y2 = safe_standardize(exp(-Y_coord % Y_coord));
    Z_exponential = arma::join_horiz(XWX_mat, arma::join_horiz(exp_x2, exp_y2));
  } else {
    Z_original = arma::join_horiz(X_coord, Y_coord);
    Z_periodic = arma::join_horiz(cos(2 * M_PI * X_coord), cos(2 * M_PI * Y_coord));

    arma::vec exp_x2 = safe_standardize(exp(-X_coord % X_coord));
    arma::vec exp_y2 = safe_standardize(exp(-Y_coord % Y_coord));
    Z_exponential = arma::join_horiz(exp_x2, exp_y2);
  }

  // Number of groups G = max(group) + 1 (group indices start at 1).
  int G = 0;
  for (int i = 0; i < p; ++i) {
    Rcpp::IntegerVector groups = gene_group_list[i];
    for (int j = 0; j < groups.length(); ++j) {
      if (groups[j] > G) {
        G = groups[j];
      }
    }
  }
  G = G + 1;

  Rcout << "Number of groups (G): " << G << std::endl;

  arma::mat WY = weights * Y;

  Rcpp::List results_original = run_mcmc(Y, n, p, gene_group_list, G, Z_original, WY,
                                         iter, burn, a_rho, b_rho, a_tau, b_tau,
                                         a_gamma, b_gamma, a_sigma, b_sigma,
                                         precomputedTraces, order, m, 0.0, "linear",
                                         informative_gamma, a_gamma_on, b_gamma_on,
                                         a_gamma_off, b_gamma_off);

  Rcpp::List results_periodic = run_mcmc(Y, n, p, gene_group_list, G, Z_periodic, WY,
                                         iter, burn, a_rho, b_rho, a_tau, b_tau,
                                         a_gamma, b_gamma, a_sigma, b_sigma,
                                         precomputedTraces, order, m, 0.0, "periodic",
                                         informative_gamma, a_gamma_on, b_gamma_on,
                                         a_gamma_off, b_gamma_off);

  Rcpp::List results_exponential = run_mcmc(Y, n, p, gene_group_list, G, Z_exponential, WY,
                                            iter, burn, a_rho, b_rho, a_tau, b_tau,
                                            a_gamma, b_gamma, a_sigma, b_sigma,
                                            precomputedTraces, order, m, 0.0, "focal",
                                            informative_gamma, a_gamma_on, b_gamma_on,
                                            a_gamma_off, b_gamma_off);

  double avg_loglike_original = results_original["avg_log_likelihood"];
  double avg_loglike_periodic = results_periodic["avg_log_likelihood"];
  double avg_loglike_exponential = results_exponential["avg_log_likelihood"];

  std::string selected_Z;
  if (avg_loglike_original >= avg_loglike_periodic && avg_loglike_original >= avg_loglike_exponential) {
    selected_Z = "original";
  } else if (avg_loglike_periodic >= avg_loglike_original && avg_loglike_periodic >= avg_loglike_exponential) {
    selected_Z = "periodic";
  } else {
    selected_Z = "exponential";
  }

  return Rcpp::List::create(
    Rcpp::Named("selected_Z") = selected_Z,
    Rcpp::Named("original_result") = results_original,
    Rcpp::Named("periodic_result") = results_periodic,
    Rcpp::Named("exponential_result") = results_exponential
  );
}
