// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
#include <omp.h>
#include <chrono>


using namespace Rcpp;
using namespace arma;
using namespace std::chrono;





double computeSingleTraceEstimate(const arma::sp_mat &W, int i, const arma::vec &v) {
  arma::vec Wv = v; 
  for (int k = 0; k < i; ++k) {
    Wv = W * Wv;
  }
  return arma::dot(v, Wv); 
}


double approximateLogDetPrecomputed(const std::vector<std::vector<double>>& precomputedTraces, double rho, int order, int m) {
  double logDetEstimate = 0.0;
  for (int i = 1; i <= order; ++i) {
    double traceEstimateSumForOrder = 0.0; 
    for (int j = 0; j < m; ++j) {
      traceEstimateSumForOrder += precomputedTraces[i-1][j]; 
    }
    logDetEstimate -= std::pow(rho, i) * traceEstimateSumForOrder / i; 
  }
  return logDetEstimate / m;
}




arma::vec log_likelihood(const arma::mat& y_norm, const arma::vec &rho, const arma::vec &sigma,const arma::mat &Z,
                         const arma::mat &psi,const arma::mat &WY,
                         const int &p,const int &n,const int &m,const int &order,const std::vector<std::vector<double>>&precomputedTraces) {
  arma::vec log_likes(p); 
  for (int j = 0; j < y_norm.n_cols; ++j) {
    double det_part = approximateLogDetPrecomputed(precomputedTraces, rho(j) , order,m);
    double sgm_part = -(n/2)*log(2*M_PI*sigma(j));
    
    arma::vec e_part = y_norm.col(j) - rho(j) * WY.col(j) - Z * psi.col(j);
    double exp_part = -arma::dot(e_part, e_part) / (2 * sigma(j));
    log_likes(j) = det_part + sgm_part + exp_part;
    
  }
  
  
  
  return log_likes;
}



// [[Rcpp::export]]
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
                    const std::string& model_name){
  

  int  j, it;
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
  arma::vec tau1_product(p, fill::zeros);  
  arma::vec tau2_product(p, fill::zeros);  
  arma::ivec gamma_1(p, fill::zeros);
  arma::ivec gamma_2(p, fill::zeros);
  
  
  arma::vec tau_1_sum(G, fill::zeros);
  arma::vec tau_2_sum(G, fill::zeros);
  arma::vec tau1_product_sum(p, fill::zeros);
  arma::vec tau2_product_sum(p, fill::zeros);
  arma::vec gamma_1_sum(p, fill::zeros);
  arma::vec gamma_2_sum(p, fill::zeros);
  
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
  

  double time_tau1_ms = 0.0;
  double time_tau2_ms = 0.0;
  double time_gamma1_ms = 0.0;
  double time_gamma2_ms = 0.0;
  double time_psi_ms = 0.0;
  double time_sigma_ms = 0.0;
  double time_rho_ms = 0.0;

  double time_compute_ZtZ_ms = 0.0;
  double time_group_prep_ms = 0.0;
  double time_tauprod_ms = 0.0;
  double time_loglike_ms = 0.0;
  double time_store_ms = 0.0;

  for (int g = 0; g < G; g++) {
    tau_1(g) = 1;  
    tau_2(g) = 1; 
  }
  

  for (int j = 0; j < p; j++) {
    Rcpp::IntegerVector gene_groups = gene_group_list[j];
    
    double tau1_effect = 1.0;
    double tau2_effect = 1.0;
    for (int k = 0; k < gene_groups.length(); ++k) {
      tau1_effect *= tau_1[gene_groups[k]];
      tau2_effect *= tau_2[gene_groups[k]];
    }
    tau1_product(j) = tau1_effect;
    tau2_product(j) = tau2_effect;
  }
  
  
  
  
  arma::mat y_norm=Y;
  
  
  
  auto t_ztz_start = high_resolution_clock::now();
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
  auto t_ztz_end = high_resolution_clock::now();
  time_compute_ZtZ_ms += duration_cast<duration<double, std::milli>>(t_ztz_end - t_ztz_start).count();
  arma::uword base_cols = Z.n_cols - 2;
  arma::mat Z_base = Z.cols(0, base_cols - 1);

  std::vector<std::vector<int>> group_genes(G);
  std::vector<std::vector<int>> gene_groups_vec(p);
  

  auto t_group_prep_start = high_resolution_clock::now();
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
  auto t_group_prep_end = high_resolution_clock::now();
  time_group_prep_ms += duration_cast<duration<double, std::milli>>(t_group_prep_end - t_group_prep_start).count();
  

  Rcout << "[Model] " << model_name ;
  int checkpoint_iter = std::max(1, iter / 5);
  
  // MCMC
  for (it = 0; it < iter; it++) {
    arma::mat tildeY_all(y_norm.n_rows, p);
    arma::mat R_base_all(y_norm.n_rows, p);
    for (int j_pre = 0; j_pre < p; ++j_pre) {
      tildeY_all.col(j_pre) = y_norm.col(j_pre) - rho(j_pre) * WY.col(j_pre);
      R_base_all.col(j_pre) = Z_base * psi.col(j_pre).rows(0, base_cols - 1);
    }
    

    auto t_tauprod_start = high_resolution_clock::now();
    for (int j = 0; j < p; ++j) {
      Rcpp::IntegerVector gene_groups = gene_group_list[j];
      

      double tau1_effect = 1.0;
      for (int k = 0; k < gene_groups.length(); ++k) {
        tau1_effect *= tau_1[gene_groups[k]];
      }
      tau1_product(j) = tau1_effect;
      

      double tau2_effect = 1.0;
      for (int k = 0; k < gene_groups.length(); ++k) {
        tau2_effect *= tau_2[gene_groups[k]];
      }
      tau2_product(j) = tau2_effect;
    }
    auto t_tauprod_end = high_resolution_clock::now();
    time_tauprod_ms += duration_cast<duration<double, std::milli>>(t_tauprod_end - t_tauprod_start).count();
    
    
    
    if (ZtZ.n_rows > 2) {  
      
      // tau_1
      auto t1_start = high_resolution_clock::now();
      std::vector<double> posterior_probs_tau1(G);
      
      for (int g_idx = 0; g_idx < G; ++g_idx) {
        double pi_k1 = double(a_tau) / (a_tau + b_tau);
        double log_prob_tau_k1_1 = log(pi_k1);
        double log_prob_tau_k1_0 = log(1 - pi_k1);
        double log_likelihood_tau_k1_1 = 0.0;
        double log_likelihood_tau_k1_0 = 0.0;
        
        const std::vector<int>& genes_in_group = group_genes[g_idx];
#pragma omp parallel for
        for (int j : genes_in_group) {
          const std::vector<int>& gene_groups = gene_groups_vec[j];
          double other_groups_effect = 1.0;
          for (const int& group : gene_groups) {
            if (group != g_idx && group >= 0 && group < G) {
              other_groups_effect *= tau_1[group];
            }
          }
          
          arma::vec tilde_y_j = tildeY_all.col(j);
          
          arma::vec e_j_1 = tilde_y_j - R_base_all.col(j)
            - (other_groups_effect * 1) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
            - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
            
            arma::vec e_j_0 = tilde_y_j - R_base_all.col(j)
              - (other_groups_effect * 0) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
              - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
              
              log_likelihood_tau_k1_1 -= arma::dot(e_j_1, e_j_1) / (2 * sigma(j));
              log_likelihood_tau_k1_0 -= arma::dot(e_j_0, e_j_0) / (2 * sigma(j));
        }
        
        log_prob_tau_k1_1 += log_likelihood_tau_k1_1;
        log_prob_tau_k1_0 += log_likelihood_tau_k1_0;
        
        double max_log_prob = std::max(log_prob_tau_k1_1, log_prob_tau_k1_0);
        double prob_tau_k1_1 = exp(log_prob_tau_k1_1 - max_log_prob);
        double prob_tau_k1_0 = exp(log_prob_tau_k1_0 - max_log_prob);
        
        posterior_probs_tau1[g_idx] = prob_tau_k1_1 / (prob_tau_k1_0 + prob_tau_k1_1);
        tau_1(g_idx) = R::rbinom(1, posterior_probs_tau1[g_idx]);
      }
      auto t1_end = high_resolution_clock::now();
      time_tau1_ms += duration_cast<duration<double, std::milli>>(t1_end - t1_start).count();

      
      // tau_2
      auto t2_start = high_resolution_clock::now();
      std::vector<double> posterior_probs_tau2(G);
      
      for (int g_idx = 0; g_idx < G; ++g_idx) {
        double pi_k2 = double(a_tau) / (a_tau + b_tau);
        double log_prob_tau_k2_1 = log(pi_k2);
        double log_prob_tau_k2_0 = log(1 - pi_k2);
        double log_likelihood_tau_k2_1 = 0.0;
        double log_likelihood_tau_k2_0 = 0.0;
        
        const std::vector<int>& genes_in_group = group_genes[g_idx];
#pragma omp parallel for
        for (int j : genes_in_group) {
          const std::vector<int>& gene_groups = gene_groups_vec[j];
          double other_groups_effect = 1.0;
          for (const int& group : gene_groups) {
            if (group != g_idx && group >= 0 && group < G) {
              other_groups_effect *= tau_2[group];
            }
          }
          
          arma::vec tilde_y_j = tildeY_all.col(j);
          
          arma::vec e_j_1 = tilde_y_j - R_base_all.col(j)
            - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
            - (other_groups_effect * 1) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
            
            arma::vec e_j_0 = tilde_y_j - R_base_all.col(j)
              - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
              - (other_groups_effect * 0) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
              
              log_likelihood_tau_k2_1 -= arma::dot(e_j_1, e_j_1) / (2 * sigma(j));
              log_likelihood_tau_k2_0 -= arma::dot(e_j_0, e_j_0) / (2 * sigma(j));
        }
        
        log_prob_tau_k2_1 += log_likelihood_tau_k2_1;
        log_prob_tau_k2_0 += log_likelihood_tau_k2_0;
        
        double max_log_prob = std::max(log_prob_tau_k2_1, log_prob_tau_k2_0);
        double prob_tau_k2_1 = exp(log_prob_tau_k2_1 - max_log_prob);
        double prob_tau_k2_0 = exp(log_prob_tau_k2_0 - max_log_prob);
        
        posterior_probs_tau2[g_idx] = prob_tau_k2_1 / (prob_tau_k2_0 + prob_tau_k2_1);
        tau_2(g_idx) = R::rbinom(1, posterior_probs_tau2[g_idx]);
      }
      auto t2_end = high_resolution_clock::now();
      time_tau2_ms += duration_cast<duration<double, std::milli>>(t2_end - t2_start).count();
      
      
      
      // gamma_1
      auto tg1_start = high_resolution_clock::now();
#pragma omp parallel for
      for (int j = 0; j < p; ++j) {
        double pi_j1 = double(a_gamma) / (a_gamma + b_gamma);  
        double log_prob_gamma_j1_1 = log(pi_j1);  
        double log_prob_gamma_j1_0 = log(1 - pi_j1);
        arma::vec tilde_y_jk = tildeY_all.col(j);
        arma::vec e_j1_1 = tilde_y_jk - R_base_all.col(j)- tau1_product(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)- tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
        arma::vec e_j1_0 = tilde_y_jk - R_base_all.col(j) - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
        log_prob_gamma_j1_1 -= arma::dot(e_j1_1, e_j1_1) / (2 * sigma(j));
        log_prob_gamma_j1_0 -= arma::dot(e_j1_0, e_j1_0) / (2 * sigma(j));
        double max_log_prob = std::max(log_prob_gamma_j1_1, log_prob_gamma_j1_0);
        double prob_gamma_j1_1 = exp(log_prob_gamma_j1_1 - max_log_prob);
        double prob_gamma_j1_0 = exp(log_prob_gamma_j1_0 - max_log_prob);
        double posterior_prob_gamma_j1_1 = prob_gamma_j1_1 / (prob_gamma_j1_1 + prob_gamma_j1_0);
        gamma_1(j) = R::rbinom(1, posterior_prob_gamma_j1_1);
      }
      auto tg1_end = high_resolution_clock::now();
      time_gamma1_ms += duration_cast<duration<double, std::milli>>(tg1_end - tg1_start).count();
      
      
      
      // gamma_2
      auto tg2_start = high_resolution_clock::now();
#pragma omp parallel for
      for (int j = 0; j < p; ++j) {
        double pi_j2 = double(a_gamma) / (a_gamma + b_gamma); 
        double log_prob_gamma_j2_1 = log(pi_j2);  
        double log_prob_gamma_j2_0 = log(1 - pi_j2);  
        arma::vec tilde_y_jk = tildeY_all.col(j);
        arma::vec e_j2_1 = tilde_y_jk - R_base_all.col(j)- tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)- tau2_product(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
        arma::vec e_j2_0 = tilde_y_jk - R_base_all.col(j)- tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j);
        log_prob_gamma_j2_1 -= arma::dot(e_j2_1, e_j2_1) / (2 * sigma(j));
        log_prob_gamma_j2_0 -= arma::dot(e_j2_0, e_j2_0) / (2 * sigma(j));
        double max_log_prob = std::max(log_prob_gamma_j2_1, log_prob_gamma_j2_0);
        double prob_gamma_j2_1 = exp(log_prob_gamma_j2_1 - max_log_prob);
        double prob_gamma_j2_0 = exp(log_prob_gamma_j2_0 - max_log_prob);
        double posterior_prob_gamma_j2_1 = prob_gamma_j2_1 / (prob_gamma_j2_1 + prob_gamma_j2_0);
        gamma_2(j) = R::rbinom(1, posterior_prob_gamma_j2_1);
      }
      auto tg2_end = high_resolution_clock::now();
      time_gamma2_ms += duration_cast<duration<double, std::milli>>(tg2_end - tg2_start).count();
      
    }else{
      
      
      // tau_1
      auto t1_start = high_resolution_clock::now();
      std::vector<double> posterior_probs_tau1(G);
      
      for (int g_idx = 0; g_idx < G; ++g_idx) {
        double pi_k1 = double(a_tau) / (a_tau + b_tau);
        double log_prob_tau_k1_1 = log(pi_k1);
        double log_prob_tau_k1_0 = log(1 - pi_k1);
        double log_likelihood_tau_k1_1 = 0.0;
        double log_likelihood_tau_k1_0 = 0.0;
        
        const std::vector<int>& genes_in_group = group_genes[g_idx];
#pragma omp parallel for
        for (int j : genes_in_group) {
          const std::vector<int>& gene_groups = gene_groups_vec[j];
          double other_groups_effect = 1.0;
          for (const int& group : gene_groups) {
            if (group != g_idx && group >= 0 && group < G) {
              other_groups_effect *= tau_1[group];
            }
          }
          
          arma::vec tilde_y_j = Y.col(j) - rho(j) * WY.col(j);
          
          arma::vec e_j_1 = tilde_y_j 
            - (other_groups_effect * 1) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
            - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
            
            arma::vec e_j_0 = tilde_y_j 
              - (other_groups_effect * 0) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
              - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
              
              log_likelihood_tau_k1_1 -= arma::dot(e_j_1, e_j_1) / (2 * sigma(j));
              log_likelihood_tau_k1_0 -= arma::dot(e_j_0, e_j_0) / (2 * sigma(j));
        }
        
        log_prob_tau_k1_1 += log_likelihood_tau_k1_1;
        log_prob_tau_k1_0 += log_likelihood_tau_k1_0;
        
        double max_log_prob = std::max(log_prob_tau_k1_1, log_prob_tau_k1_0);
        double prob_tau_k1_1 = exp(log_prob_tau_k1_1 - max_log_prob);
        double prob_tau_k1_0 = exp(log_prob_tau_k1_0 - max_log_prob);
        
        posterior_probs_tau1[g_idx] = prob_tau_k1_1 / (prob_tau_k1_0 + prob_tau_k1_1);
        tau_1(g_idx) = R::rbinom(1, posterior_probs_tau1[g_idx]);
      }
      auto t1_end = high_resolution_clock::now();
      time_tau1_ms += duration_cast<duration<double, std::milli>>(t1_end - t1_start).count();

      // tau_2
      auto t2_start = high_resolution_clock::now();
      std::vector<double> posterior_probs_tau2(G);
      
      for (int g_idx = 0; g_idx < G; ++g_idx) {
        double pi_k2 = double(a_tau) / (a_tau + b_tau);
        double log_prob_tau_k2_1 = log(pi_k2);
        double log_prob_tau_k2_0 = log(1 - pi_k2);
        double log_likelihood_tau_k2_1 = 0.0;
        double log_likelihood_tau_k2_0 = 0.0;
        
        const std::vector<int>& genes_in_group = group_genes[g_idx];
#pragma omp parallel for
        for (int j : genes_in_group) {
          const std::vector<int>& gene_groups = gene_groups_vec[j];
          double other_groups_effect = 1.0;
          for (const int& group : gene_groups) {
            if (group != g_idx && group >= 0 && group < G) {
              other_groups_effect *= tau_2[group];
            }
          }
          
          arma::vec tilde_y_j = Y.col(j) - rho(j) * WY.col(j);
          
          arma::vec e_j_1 = tilde_y_j 
            - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
            - (other_groups_effect * 1) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
            
            arma::vec e_j_0 = tilde_y_j 
              - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)
              - (other_groups_effect * 0) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
              
              log_likelihood_tau_k2_1 -= arma::dot(e_j_1, e_j_1) / (2 * sigma(j));
              log_likelihood_tau_k2_0 -= arma::dot(e_j_0, e_j_0) / (2 * sigma(j));
        }
        
        log_prob_tau_k2_1 += log_likelihood_tau_k2_1;
        log_prob_tau_k2_0 += log_likelihood_tau_k2_0;
        
        double max_log_prob = std::max(log_prob_tau_k2_1, log_prob_tau_k2_0);
        double prob_tau_k2_1 = exp(log_prob_tau_k2_1 - max_log_prob);
        double prob_tau_k2_0 = exp(log_prob_tau_k2_0 - max_log_prob);
        
        posterior_probs_tau2[g_idx] = prob_tau_k2_1 / (prob_tau_k2_0 + prob_tau_k2_1);
        tau_2(g_idx) = R::rbinom(1, posterior_probs_tau2[g_idx]);
      }
      auto t2_end = high_resolution_clock::now();
      time_tau2_ms += duration_cast<duration<double, std::milli>>(t2_end - t2_start).count();
      

      
      // gamma_1
      auto tg1_start = high_resolution_clock::now();
#pragma omp parallel for
      for (int j = 0; j < p; ++j) {
        double pi_j1 = double(a_gamma) / (a_gamma + b_gamma);  
        double log_prob_gamma_j1_1 = log(pi_j1);  
        double log_prob_gamma_j1_0 = log(1 - pi_j1);
        arma::vec tilde_y_jk = Y.col(j) - rho(j) * WY.col(j);
        arma::vec e_j1_1 = tilde_y_jk - tau1_product(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)- tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
        arma::vec e_j1_0 = tilde_y_jk  - tau2_product(j) * gamma_2(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
        log_prob_gamma_j1_1 -= arma::dot(e_j1_1, e_j1_1) / (2 * sigma(j));
        log_prob_gamma_j1_0 -= arma::dot(e_j1_0, e_j1_0) / (2 * sigma(j));
        double max_log_prob = std::max(log_prob_gamma_j1_1, log_prob_gamma_j1_0);
        double prob_gamma_j1_1 = exp(log_prob_gamma_j1_1 - max_log_prob);
        double prob_gamma_j1_0 = exp(log_prob_gamma_j1_0 - max_log_prob);
        double posterior_prob_gamma_j1_1 = prob_gamma_j1_1 / (prob_gamma_j1_1 + prob_gamma_j1_0);
        gamma_1(j) = R::rbinom(1, posterior_prob_gamma_j1_1);
      }
      auto tg1_end = high_resolution_clock::now();
      time_gamma1_ms += duration_cast<duration<double, std::milli>>(tg1_end - tg1_start).count();
      
      
      
      // gamma_2
      auto tg2_start = high_resolution_clock::now();
#pragma omp parallel for
      for (int j = 0; j < p; ++j) {
        double pi_j2 = double(a_gamma) / (a_gamma + b_gamma); 
        double log_prob_gamma_j2_1 = log(pi_j2);  
        double log_prob_gamma_j2_0 = log(1 - pi_j2);  
        arma::vec tilde_y_jk = Y.col(j) - rho(j) * WY.col(j);
        arma::vec e_j2_1 = tilde_y_jk - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j)- tau2_product(j) * Z.col(Z.n_cols - 1) * psi(psi.n_rows - 1, j);
        arma::vec e_j2_0 = tilde_y_jk - tau1_product(j) * gamma_1(j) * Z.col(Z.n_cols - 2) * psi(psi.n_rows - 2, j);
        log_prob_gamma_j2_1 -= arma::dot(e_j2_1, e_j2_1) / (2 * sigma(j));
        log_prob_gamma_j2_0 -= arma::dot(e_j2_0, e_j2_0) / (2 * sigma(j));
        double max_log_prob = std::max(log_prob_gamma_j2_1, log_prob_gamma_j2_0);
        double prob_gamma_j2_1 = exp(log_prob_gamma_j2_1 - max_log_prob);
        double prob_gamma_j2_0 = exp(log_prob_gamma_j2_0 - max_log_prob);
        double posterior_prob_gamma_j2_1 = prob_gamma_j2_1 / (prob_gamma_j2_1 + prob_gamma_j2_0);
        gamma_2(j) = R::rbinom(1, posterior_prob_gamma_j2_1);
      }
      auto tg2_end = high_resolution_clock::now();
      time_gamma2_ms += duration_cast<duration<double, std::milli>>(tg2_end - tg2_start).count();
      
      
      
      
    }
    
    
    
    // psi
    auto tpsi_start = high_resolution_clock::now();
#pragma omp parallel for
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
    }
    auto tpsi_end = high_resolution_clock::now();
    time_psi_ms += duration_cast<duration<double, std::milli>>(tpsi_end - tpsi_start).count();
    
    
    //sigma
    auto tsigma_start = high_resolution_clock::now();
#pragma omp parallel for
    for (int j = 0; j < p; j++) {
      arma::vec Z_psi_part = Z * psi.col(j);
      Z_psi_part(Z.n_cols - 2) *= tau1_product(j) * gamma_1(j);
      Z_psi_part(Z.n_cols - 1) *= tau2_product(j) * gamma_2(j);
      arma::vec e_kj = y_norm.col(j) - rho(j) * WY.col(j)  - Z_psi_part;
      if (e_kj.has_nan() || e_kj.has_inf()) {
        Rcout << "e_kj has invalid values for j: " << j << std::endl;
        continue;
      }
      double alpha = n / 2.0 + a_sigma / 2.0;
      double beta = arma::dot(e_kj, e_kj) / 2.0 + b_sigma / 2.0;
      sigma(j) = 1 / R::rgamma(alpha, 1/beta);
    }
    auto tsigma_end = high_resolution_clock::now();
    time_sigma_ms += duration_cast<duration<double, std::milli>>(tsigma_end - tsigma_start).count();
    
    //rho
    auto trho_start = high_resolution_clock::now();
#pragma omp parallel for
    for (int j = 0; j < p; j++) {
      hastings = 0;
      double rho_temp = rho(j)+ c[j] * R::rnorm(0, 1);
      double det_current = approximateLogDetPrecomputed(precomputedTraces, rho(j) , order,m);
      double det_proposed = approximateLogDetPrecomputed(precomputedTraces, rho_temp , order,m);
      double det_ratio = exp(det_proposed-det_current);
      arma::vec psi_corrected = psi.col(j);
      arma::vec e_kj_old = y_norm.col(j) - Z * psi.col(j) - rho(j) * WY.col(j);
      arma::vec e_kj_new = y_norm.col(j) - Z * psi.col(j) - rho_temp * WY.col(j);
      double error_ratio = exp(-arma::dot(e_kj_new, e_kj_new) / (2 * sigma(j)) + arma::dot(e_kj_old, e_kj_old) / (2 * sigma(j)));
      hastings += error_ratio * det_ratio;
      if (hastings >= double(rand() % 10001) / 10000 && rho_temp > -1 && rho_temp < 1) {
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
    auto trho_end = high_resolution_clock::now();
    time_rho_ms += duration_cast<duration<double, std::milli>>(trho_end - trho_start).count();
    
    

    if ((it+1) % checkpoint_iter == 0 || (it+1) == iter) {
      int pct = (it+1) * 100 / iter;
      Rcout << "[Model] " << model_name << " progress: " << pct << "%\n";
    }

    if (it >= burn) {
      int store_index = it - burn;
   
      for (j = 0; j < p; j++) {
        
         rho_sum(j) += rho(j);
         sigma_sum(j) += sigma(j);
         psi_sum[j] += psi.col(j);
         
      
        gamma_1_sum(j) += gamma_1(j);
        gamma_2_sum(j) += gamma_2(j);
        
        tau1_product_sum(j) += tau1_product(j);
        tau2_product_sum(j) += tau2_product(j);
        
        tau_gamma_1_sum(j) += tau1_product(j) * gamma_1(j);
        tau_gamma_2_sum(j) += tau2_product(j) * gamma_2(j);
      }
      auto t_loglike_start = high_resolution_clock::now();
      arma::vec gene_loglike = log_likelihood(y_norm, rho, sigma, Z, psi, WY, p, n, m, order, precomputedTraces);
      auto t_loglike_end = high_resolution_clock::now();
      time_loglike_ms += duration_cast<duration<double, std::milli>>(t_loglike_end - t_loglike_start).count();

      auto t_store_start = high_resolution_clock::now();
      for(int j = 0; j < p; j++) {
        loglike_sum(j) = loglike_sum(j) + gene_loglike(j);
      }
      auto t_store_end = high_resolution_clock::now();
      time_store_ms += duration_cast<duration<double, std::milli>>(t_store_end - t_store_start).count();
    }
  }
  
  
  
  
   arma::vec rho_mean = rho_sum / static_cast<double>(effective_iterations);
   arma::vec sigma_mean = sigma_sum / static_cast<double>(effective_iterations);
   for (int j = 0; j < p; j++) {
   psi_sum[j] /= static_cast<double>(effective_iterations);
   }
   
  
  
  arma::vec tau_gamma_1_mean = tau_gamma_1_sum / static_cast<double>(effective_iterations);
  arma::vec tau_gamma_2_mean = tau_gamma_2_sum / static_cast<double>(effective_iterations);
  arma::vec tau_1_mean = tau1_product_sum / static_cast<double>(effective_iterations);
  arma::vec tau_2_mean = tau2_product_sum / static_cast<double>(effective_iterations);
  arma::vec gamma_1_mean = gamma_1_sum / static_cast<double>(effective_iterations);
  arma::vec gamma_2_mean = gamma_2_sum / static_cast<double>(effective_iterations);
  
 
  
  
  
  arma::vec tau_gamma(p, arma::fill::zeros);
  for (int i = 0; i < p; ++i) {
    tau_gamma(i) = std::max(tau_gamma_1_mean(i), tau_gamma_2_mean(i));
  }
  
  
  
  
  double total_loglike_sum = 0.0;                     
  double total_loglike_mean = 0.0;                     
  
  for(int j = 0; j < p; j++) {
    loglike_gene_mean(j) = loglike_sum(j) / static_cast<double>(effective_iterations);
    total_loglike_sum += loglike_sum(j);  
  }
  

  total_loglike_mean = total_loglike_sum / static_cast<double>(effective_iterations);
  
  Rcpp::List timings_ms = Rcpp::List::create(
    Rcpp::Named("tau1") = time_tau1_ms,
    Rcpp::Named("tau2") = time_tau2_ms,
    Rcpp::Named("gamma1") = time_gamma1_ms,
    Rcpp::Named("gamma2") = time_gamma2_ms,
    Rcpp::Named("sigma") = time_sigma_ms,
    Rcpp::Named("rho") = time_rho_ms,
    Rcpp::Named("build_WY") = time_build_WY_ms,
    Rcpp::Named("compute_ZtZ") = time_compute_ZtZ_ms,
    Rcpp::Named("group_prep") = time_group_prep_ms,
    Rcpp::Named("tau_product") = time_tauprod_ms,
    Rcpp::Named("loglike") = time_loglike_ms,
    Rcpp::Named("store") = time_store_ms
  );
  Rcpp::List timings_ms_per_iter = Rcpp::List::create(
    Rcpp::Named("tau1") = time_tau1_ms / iter,
    Rcpp::Named("tau2") = time_tau2_ms / iter,
    Rcpp::Named("gamma1") = time_gamma1_ms / iter,
    Rcpp::Named("gamma2") = time_gamma2_ms / iter,
    Rcpp::Named("psi") = time_psi_ms / iter,
    Rcpp::Named("sigma") = time_sigma_ms / iter,
    Rcpp::Named("rho") = time_rho_ms / iter,
    Rcpp::Named("build_WY") = time_build_WY_ms / iter,
    Rcpp::Named("compute_ZtZ") = time_compute_ZtZ_ms / iter,
    Rcpp::Named("group_prep") = time_group_prep_ms / iter,
    Rcpp::Named("tau_product") = time_tauprod_ms / iter,
    Rcpp::Named("loglike") = time_loglike_ms / iter,
    Rcpp::Named("store") = time_store_ms / iter
  );
  
  return Rcpp::List::create(

    Rcpp::Named("tau_gamma_upper") = tau_gamma,
    Rcpp::Named("gene_log_likelihood_mean") = loglike_gene_mean,
    Rcpp::Named("avg_log_likelihood") = total_loglike_mean,
    

    Rcpp::Named("rho_mean") = rho_mean,
    Rcpp::Named("sigma_mean") = sigma_mean,
    Rcpp::Named("psi_mean") = psi_sum
    

    );
}







// [[Rcpp::export]]
Rcpp::List Baydubig_mcmc(const arma::mat& Y, 
                         const Rcpp::List& gene_group_list,  
                         const Rcpp::Nullable<arma::mat>& XWX,
                         const arma::vec& X_coord, 
                         const arma::vec& Y_coord,
                         const arma::sp_mat& weights, 
                         int iter, 
                         int burn,
                         double a_rho=1.00, 
                         double b_rho=1.00, 
                         double a_tau=1, 
                         double b_tau=1, 
                         double a_gamma=1, 
                         double b_gamma=1,
                         double a_sigma=0.01, 
                         double b_sigma=0.01) {

  int max_threads = omp_get_max_threads();
  int num_threads_to_use = max_threads;
  if (num_threads_to_use < 1) {
    num_threads_to_use = 1;
  }
  
  omp_set_num_threads(num_threads_to_use);
  
  Rcout << "OpenMP is enabled. Number of threads: " << omp_get_max_threads() << std::endl;
  Rcout << "Using " << num_threads_to_use << " threads for parallel execution." << std::endl;
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  
  int order = 5;
  int m = 1;
  
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
    // If XWX is provided, dereference it and join it with coordinates and transformations
    arma::mat XWX_mat = Rcpp::as<arma::mat>(XWX);  // Convert from Nullable<arma::mat> to arma::mat
    Z_periodic = arma::join_horiz(XWX_mat, arma::join_horiz(cos(2 * M_PI * X_coord), cos(2 * M_PI * Y_coord)));
    Z_original = arma::join_horiz(XWX_mat, arma::join_horiz(X_coord, Y_coord));
    
    Z_exponential = arma::join_horiz(XWX_mat, arma::join_horiz(exp(-X_coord % X_coord), exp(-Y_coord % Y_coord)));
  } else {
    // If XWX is not provided, just use coordinates and transformations
    Z_periodic = arma::join_horiz(cos(2 * M_PI * X_coord), cos(2 * M_PI * Y_coord));
    Z_original = arma::join_horiz(X_coord, Y_coord);
    
    Z_exponential = arma::join_horiz(exp(-X_coord % X_coord), exp(-Y_coord % Y_coord));
  }
  

  int G = 0;

  for (int i = 0; i < p; ++i) {
    Rcpp::IntegerVector groups = gene_group_list[i];
    for (int j = 0; j < groups.length(); ++j) {
      if (groups[j] > G) {
        G = groups[j];
      }
    }
  }
  
  std::cout << "Maximum group number (G): " << G << std::endl;

  Rcpp::List results_original, results_periodic, results_exponential;
  
  

  auto t_wy_start = high_resolution_clock::now();
  arma::mat WY = weights * Y;
  auto t_wy_end = high_resolution_clock::now();
  double time_build_WY_ms = duration_cast<duration<double, std::milli>>(t_wy_end - t_wy_start).count();

  results_original = run_mcmc(Y, n, p, gene_group_list, G, Z_original, WY,
                              iter, burn, a_rho, b_rho, a_tau, b_tau, a_gamma, b_gamma, 
                              a_sigma, b_sigma, precomputedTraces, order, m, time_build_WY_ms, "linear");
  
  results_periodic = run_mcmc(Y, n, p, gene_group_list, G, Z_periodic, WY,
                              iter, burn, a_rho, b_rho, a_tau, b_tau, a_gamma, b_gamma,
                              a_sigma, b_sigma, precomputedTraces, order, m, time_build_WY_ms, "periodic");
  
  results_exponential = run_mcmc(Y, n, p, gene_group_list, G, Z_exponential, WY,
                                 iter, burn, a_rho, b_rho, a_tau, b_tau, a_gamma, b_gamma,
                                 a_sigma, b_sigma, precomputedTraces, order, m, time_build_WY_ms, "focal");
  
  
  
  double avg_loglike_original = results_original["avg_log_likelihood"];
  double avg_loglike_periodic = results_periodic["avg_log_likelihood"];
  double avg_loglike_exponential = results_exponential["avg_log_likelihood"];
  
  std::string selected_Z;
  arma::mat selected_Z_mat;
  
  if (avg_loglike_original >= avg_loglike_periodic && avg_loglike_original >= avg_loglike_exponential) {
    selected_Z = "original";
    selected_Z_mat = Z_original;
  } else if (avg_loglike_periodic >= avg_loglike_original && avg_loglike_periodic >= avg_loglike_exponential) {
    selected_Z = "periodic";
    selected_Z_mat = Z_periodic;
  } else {
    selected_Z = "exponential";
    selected_Z_mat = Z_exponential;
  }
  
 
  return Rcpp::List::create(

    Rcpp::Named("selected_Z") = selected_Z,
    Rcpp::Named("original_result") = results_original,
    Rcpp::Named("periodic_result") = results_periodic,
    Rcpp::Named("exponential_result") = results_exponential
  
  );
}



