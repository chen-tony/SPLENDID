/////////////////////////////////////////////////////////////////////////////////////
/// adapted from original bigstatsr package: https://github.com/privefl/bigstatsr ///
/////////////////////////////////////////////////////////////////////////////////////

#ifndef BIGSTATSR_BIGLASSO_LIN_HPP_INCLUDED
#define BIGSTATSR_BIGLASSO_LIN_HPP_INCLUDED

#include <RcppArmadillo.h>

using std::size_t;

/******************************************************************************/

namespace slingshot { namespace slingshot_linear {

#include <utils.hpp>

using namespace slingshot::slingshot_utils;


// Gaussian loss
double gLoss(const arma::vec& r) {
  return arma::dot(r, r);
}

/******************************************************************************/
// Coordinate descent for gaussian models L0L1 (with tuning data)
template <class C>
List cdfit_gaussian_par(C macc,
                             const arma::vec& y,
                             const arma::mat& Z,
                             const arma::mat& E,
                             C macc_tun,
                             const arma::vec& y_tun,
                             const arma::mat& Z_tun,
                             const arma::mat& E_tun,
                             const arma::mat& lambda,
                             const arma::vec& Z2,
                             const arma::vec& X_scale,
                             const arma::vec& XE_scale,
                             const arma::vec& XY_norm,
                             const arma::vec& XEY_norm,
                             const Rcpp::LogicalVector& grouped,
                             const std::string & lambda_ix,
                             const std::string & pval,
                             double tol,
                             int maxit,
                             int dfmax,
                             int n_abort,
                             int nlam_min,
                             int ncheck) {
  
  size_t n = macc.nrow(); // number of observations used for fitting model
  size_t p = macc.ncol();
  size_t n_tun = macc_tun.nrow();
  int M = Z.n_cols;
  int K = E.n_cols;
  std::string message="";
  
  int L = lambda.n_rows;
  
  double metric = 0.0, metric_min = 100.0;
  int no_change = 0;
  
  arma::sp_mat out_beta(p, L*(K+1));
  arma::mat beta(p, K+1);
  // sp_mat::const_iterator beta_it;
  
  arma::rowvec beta_tilde(K+1);
  arma::rowvec beta_j(K+1);
  arma::rowvec dbeta(K+1);
  double beta_norm = 0.0;
  double L1 = 0.0;
  double L2 =0.0;
  int L0 = 0;
  int L0_het = 0;
  
  arma::vec alpha(M);
  arma::mat out_alpha(M, L);
  double alpha_m = 0.0;
  double dalpha = 0.0;
  
  arma::vec pred_tun(n_tun);
  
  // Objects to be returned to R
  Rcpp::IntegerVector iter(L, NA_INTEGER);
  Rcpp::NumericVector loss(L, NA_REAL);
  Rcpp::NumericVector metrics(L, NA_REAL);
  Rcpp::IntegerVector nb_active(L, NA_INTEGER);
  Rcpp::IntegerVector nb_interact(L, NA_INTEGER);
  Rcpp::IntegerVector nb_checks(L, NA_INTEGER);
  Rcpp::IntegerVector nb_candidate(L, NA_INTEGER);
  arma::mat pred(n_tun, L);
  
  double lambda0=0.0, lambda1=0.0, lambda2=0.0, cutoff=0.0;
  double max_update=0.0, thresh=0.0; // , shift_scaled, cpsum;
  size_t i=0, j=0, k=0, m=0, violations=0;
  double max_dbeta = 0.0;
  int check=0;
  double z=0.0; 
  
  Rcpp::LogicalVector in_A(p); // ever-active set
  Rcpp::LogicalVector in_S(p); // strong set
  int n_strong=0;
  
  // compute initial covariate effects
  alpha = arma::solve(Z, y);
  arma::vec r = y - Z * alpha;
  loss[0] = gLoss(r) / n;
  // thresh = tol * loss[0] / n;
  thresh = tol;
  
  double score0=0.0, scale0=0.0;
  arma::vec score(p); // score for strong set
  for (j=0; j<p; j++) {
    score0 = XY_norm[j] * (1-grouped[j]) + XEY_norm[j] * grouped[j];
    scale0 = X_scale[j] * (1-grouped[j]) + XE_scale[j] * grouped[j];
    if (score0 * scale0 > 0) {
      score[j] = score0 / scale0;
    } else {
      score[j] = 0;
    }
  }
  
  // Path
  for (int l = 0; l < L; l++) {
    check = 0;
    
    lambda0 = lambda(l,0);
    lambda1 = lambda(l,1);
    lambda2 = lambda(l,2);
    
    n_strong = 0;
    for (j=0; j<p; j++) {
      if (grouped[j]) {
        cutoff = std::sqrt(lambda0 / (XE_scale[j] + 2*lambda2)) + 0.5*lambda1 / (XE_scale[j] + 2*lambda2);
      } else {
        cutoff = std::sqrt(lambda0 / (X_scale[j] + 2*lambda2)) + 0.5*lambda1 / (X_scale[j] + 2*lambda2);
      }
      
      if (score[j] > cutoff) {
        in_S[j] = true;
        in_A[j] = true;
        n_strong++;
      } else {
        in_S[j] = false;
        in_A[j] = false;
      } 
    }
    
    // Approx: no check of rest set
    iter[l] = 0;
    while (iter[l] < maxit) {
      
      while (iter[l] < maxit) {
        iter[l]++;
        max_update = 0;
        
        // update covariates
        for (m=0; m<M; m++) {
          alpha_m = alpha[m]; // old value
          
          // correlation with residuals
          alpha[m] = alpha_m + arma::dot(Z.col(m), r) / Z2[m];
          
          dalpha = alpha[m] - alpha_m;
          if (std::fabs(dalpha) > 1e-20) {
            max_update = std::max(max_update, std::fabs(dalpha));
            
            // update residuals
            r -= Z.col(m) * dalpha;
          }
        }
        
        max_update = std::max(max_update, std::fabs(dalpha));
        
        // Solve lasso over ever-active set
        max_update = 0;
        L0 = 0;
        L0_het = 0;
        L1 = 0.0;
        L2 = 0.0;
        for (j = 0; j < p; j++) {
          if (in_A[j]) {
            if (grouped[j]) {
              cutoff = std::sqrt(lambda0 / (XE_scale[j] + 2*lambda2)) + 0.5*lambda1 / (XE_scale[j] + 2*lambda2);
              
              for (k=0; k<K+1; k++) {
                beta_j[k] = beta(j,k);
                beta_tilde[k] = beta(j,k);
              }
              
              for (i = 0; i < n; i++) {
                z = macc(i,j) * r[i] / n / (XE_scale[j] + 2*lambda2);
                beta_tilde[0] += z;
                // beta_tilde[0] += macc(i,j) * r[i] / n / (XE_scale[j] + 2*lambda2);
                
                for (k=1; k<K+1; k++) {
                  beta_tilde[k] += z * E(i,k-1);
                  // beta_tilde[k] += macc(i,j) * E(i,k-1) * r[i] / n / (XE_scale[j] + 2*lambda2);
                }
              }
              
              beta_norm = my_vecnorm(beta_tilde);
              score[j] = beta_norm;
              
              // thresholding
              if (beta_norm > cutoff) {
                // update beta
                beta.row(j) = beta_tilde * (1 - 0.5*lambda1 / beta_norm / (XE_scale[j] + 2*lambda2));
                
                L0 ++;
                L0_het ++; 
                L1 += beta_norm - 0.5*lambda1 / (XE_scale[j] + 2*lambda2);
                L2 += (beta_norm - 0.5*lambda1 / (XE_scale[j] + 2*lambda2)) * (beta_norm - 0.5*lambda1 / (XE_scale[j] + 2*lambda2));
                
                // difference in beta
                dbeta[0] = beta(j,0) - beta_j[0];
                max_dbeta = std::fabs(dbeta[0]);
                for (k=1; k<K+1; k++) {
                  dbeta[k] = beta(j,k) - beta_j[k];
                  max_dbeta = std::max(max_dbeta, std::fabs(dbeta[k]));
                }
              } else {
                for (k=0; k<K+1; k++) {
                  beta(j,k) = 0.0;
                  dbeta[k] = -beta_j[k];
                  max_dbeta = std::max(max_dbeta, std::fabs(beta_j[k]));
                }
                in_A[j] = false;
              }
              
            } else {
              cutoff = std::sqrt(lambda0 / (X_scale[j] + 2*lambda2)) + 0.5*lambda1 / (X_scale[j] + 2*lambda2);
              
              beta_j[0] = beta(j,0);
              beta_tilde[0] = beta(j,0);
              
              for (i = 0; i < n; i++) {
                beta_tilde[0] += macc(i,j) * r[i] / n / X_scale[j];
              }
              
              beta_norm = std::fabs(beta_tilde[0]);
              score[j] = beta_norm;
              
              // thresholding
              if (beta_norm > cutoff) {
                // update beta
                beta(j,0) = beta_tilde[0] * (1 - 0.5*lambda1 / beta_norm / (X_scale[j] + 2*lambda2));
                
                L0 ++;
                L1 += std::fabs(beta(j,0));
                L2 += beta(j,0) * beta(j,0);
                
                // difference in beta
                dbeta[0] = beta(j,0) - beta_j[0];
                max_dbeta = std::fabs(dbeta[0]);
              } else {
                beta(j,0) = 0.0;
                dbeta[0] = -beta_j[0];
                max_dbeta = std::fabs(beta_j[0]);
                in_A[j] = false;
              }
            }
            
            max_update = std::max(max_update, max_dbeta);
            
            // update residuals
            if (std::fabs(dbeta[0]) > 1e-20) {
              for (i=0; i<n; i++) {
                r(i) -= macc(i,j) * dbeta[0];
              }
            }
            if (grouped[j]) {
              for (k=1; k<K+1; k++) {
                if (std::fabs(dbeta[k]) > 1e-20) {
                  for (i=0; i<n; i++) {
                    r(i) -= macc(i,j) * E(i,k-1) * dbeta[k];
                  }
                }
              }
            }
          }
        }
        // Check for convergence
        if (max_update < thresh) {
          break;
        }
      }
      
      // Scan for violations in strong set
      violations = check_strong_set(
        in_A, in_S, macc, E, score, X_scale, XE_scale, grouped, beta, lambda0, lambda1, lambda2, r);
      if (violations == 0 || check > ncheck) {
        break;
      } else {
        check++;
      }
    }
    
    // objective
    loss[l] = gLoss(r) / n + lambda0*L0 + lambda1*L1 + lambda2*L2;
    nb_active[l] = L0;
    nb_interact[l] = L0_het;
    nb_checks[l] = check;
    nb_candidate[l] = n_strong;
    
    // prediction
    pred_tun = predict(macc_tun, beta, alpha, Z_tun, E_tun, grouped, in_A);
    metric = gLoss(pred_tun - y_tun) / n;
    // metric = arma::cor(pred_tun, y_tun);
    metrics[l] = metric;
    if (metric < metric_min) {
      metric_min = metric;
      no_change = 0;
    }
    for (i=0; i<n_tun; i++) pred(i,l) = pred_tun[i];
    
    // update output
    //// SNP effects
    // for (beta_it=beta.begin(); beta_it!=beta.end(); beta_it++) {
    //   j=beta_it.row();
    //   k=beta_it.col();
    //   
    //   out_beta(j, (K+1)*l + k) = *beta_it;
    // }
    
    for (j=0; j<p; j++) {
      if (in_A[j]) {
        for (k=0; k<K+1; k++) {
          out_beta(j, (K+1)*l + k) = beta(j,k);
        }
      }
    }
    
    //// covariate effects
    // for (m=0; m<M; m++) {
    //   out(p*(K+1) + m, l) = alpha[m];
    // }
    for (m=0; m<M; m++) {
      out_alpha(m, l) = alpha[m];
    }
    
    // check dfmax
    if (nb_active[l] >= dfmax) {
      return List::create(_["beta"] = out_beta, 
                          _["alpha"] = out_alpha, 
                          _["loss"] = loss, 
                          _["iter"] = iter, 
                          _["metrics"] = metrics, 
                          _["pred"] = pred, 
                          _["messages"] = "Too many variables",
                          _["nb_active"] = nb_active, 
                          _["nb_interact"] = nb_interact, 
                          _["nb_checks"] = nb_checks, 
                          _["nb_candidate"] = nb_candidate);
    }
    
    message = lambda_ix + "." + pval + "." + std::to_string(l) + ":" + std::to_string(iter[l]) + "," + std::to_string(L0) + "," + std::to_string(metric) + "\n";
    Rcout << message;
    if (l > 0 && metric >= metrics[l - 1]) no_change++;
    
    // check metric
    if (metric >= 10) {
      return List::create(_["beta"] = out_beta, 
                          _["alpha"] = out_alpha, 
                          _["loss"] = loss, 
                          _["iter"] = iter, 
                          _["metrics"] = metrics, 
                          //_["pred"] = pred, 
                          _["messages"] = "Diverged",
                          _["nb_active"] = nb_active, 
                          _["nb_interact"] = nb_interact, 
                          _["nb_checks"] = nb_checks, 
                          _["nb_candidate"] = nb_candidate);
    }
    
    // make submatrices // try sparse matrices
    if (l >= nlam_min && no_change >= n_abort) {
      return List::create(_["beta"] = out_beta, 
                          _["alpha"] = out_alpha, 
                          _["loss"] = loss, 
                          _["iter"] = iter, 
                          _["metrics"] = metrics, 
                          //_["pred"] = pred, 
                          _["messages"] = "No more improvement",
                          _["nb_active"] = nb_active, 
                          _["nb_interact"] = nb_interact, 
                          _["nb_checks"] = nb_checks, 
                          _["nb_candidate"] = nb_candidate);
    }
  }
  
  return List::create(_["beta"] = out_beta, 
                      _["alpha"] = out_alpha, 
                      _["loss"] = loss, 
                      _["iter"] = iter, 
                      _["metrics"] = metrics, 
                      //_["pred"] = pred, 
                      _["messages"] = "Complete path",
                      _["nb_active"] = nb_active, 
                      _["nb_interact"] = nb_interact, 
                      _["nb_checks"] = nb_checks, 
                      _["nb_candidate"] = nb_candidate);
}

/******************************************************************************/

// Coordinate descent for gaussian models L0L1 (no tuning data)
// for parallelized, minimize printing
template <class C>
List cdfit_gaussian_par0(C macc,
                              const arma::vec& y,
                              const arma::mat& Z,
                              const arma::mat& E,
                              const arma::mat& lambda,
                              const arma::vec& Z2,
                              const arma::vec& X_scale,
                              const arma::vec& XE_scale,
                              const arma::vec& XY_norm,
                              const arma::vec& XEY_norm,
                              const Rcpp::LogicalVector& grouped,
                              const std::string & lambda_ix,
                              const std::string & pval,
                              double tol,
                              int maxit,
                              int dfmax,
                              int ncheck) {
  
  size_t n = macc.nrow(); // number of observations used for fitting model
  size_t p = macc.ncol();
  int M = Z.n_cols;
  int K = E.n_cols;
  std::string message="";
  
  int L = lambda.n_rows;
  
  arma::sp_mat out_beta(p, L*(K+1));
  arma::mat beta(p, K+1);
  // sp_mat::const_iterator beta_it;
  
  arma::rowvec beta_tilde(K+1);
  arma::rowvec beta_j(K+1);
  arma::rowvec dbeta(K+1);
  double beta_norm = 0.0;
  double L1 = 0.0;
  double L2 =0.0;
  int L0 = 0;
  int L0_het = 0;
  
  arma::vec alpha(M);
  arma::mat out_alpha(M, L);
  double alpha_m = 0.0;
  double dalpha = 0.0;
  
  // Objects to be returned to R
  Rcpp::IntegerVector iter(L, NA_INTEGER);
  Rcpp::NumericVector loss(L, NA_REAL);
  Rcpp::NumericVector metrics(L, NA_REAL);
  Rcpp::IntegerVector nb_active(L, NA_INTEGER);
  Rcpp::IntegerVector nb_interact(L, NA_INTEGER);
  Rcpp::IntegerVector nb_checks(L, NA_INTEGER);
  Rcpp::IntegerVector nb_candidate(L, NA_INTEGER);
  
  double lambda0=0.0, lambda1=0.0, lambda2=0.0, cutoff=0.0;
  double max_update=0.0, thresh=0.0; // , shift_scaled, cpsum;
  size_t i=0, j=0, k=0, m=0, violations=0;
  double max_dbeta = 0.0;
  int check=0;
  double z=0.0; 
  
  Rcpp::LogicalVector in_A(p); // ever-active set
  Rcpp::LogicalVector in_S(p); // strong set
  int n_strong=0;
  
  // compute initial covariate effects
  alpha = arma::solve(Z, y);
  arma::vec r = y - Z * alpha;
  loss[0] = gLoss(r) / n;
  // thresh = tol * loss[0] / n;
  thresh = tol;
  
  double score0=0.0, scale0=0.0;
  arma::vec score(p); // score for strong set
  for (j=0; j<p; j++) {
    score0 = XY_norm[j] * (1-grouped[j]) + XEY_norm[j] * grouped[j];
    scale0 = X_scale[j] * (1-grouped[j]) + XE_scale[j] * grouped[j];
    if (score0 * scale0 > 0) {
      score[j] = score0 / scale0;
    } else {
      score[j] = 0;
    }
  }
  
  // Path
  for (int l = 0; l < L; l++) {
    
    check = 0;
    
    lambda0 = lambda(l,0);
    lambda1 = lambda(l,1);
    lambda2 = lambda(l,2);
    
    n_strong = 0;
    for (j=0; j<p; j++) {
      if (grouped[j]) {
        cutoff = std::sqrt(lambda0 / (XE_scale[j] + 2*lambda2)) + 0.5*lambda1 / (XE_scale[j] + 2*lambda2);
      } else {
        cutoff = std::sqrt(lambda0 / (X_scale[j] + 2*lambda2)) + 0.5*lambda1 / (X_scale[j] + 2*lambda2);
      }
      
      if (score[j] > cutoff) {
        in_S[j] = true;
        in_A[j] = true;
        n_strong++;
      } else {
        in_S[j] = false;
        in_A[j] = false;
        for (k=0; k<K+1; k++) beta(j,k) = 0;
      } 
    }
    
    // Approx: no check of rest set
    iter[l] = 0;
    while (iter[l] < maxit) {
      while (iter[l] < maxit) {
        
        iter[l]++;
        max_update = 0;
        
        // update covariates
        for (m=0; m<M; m++) {
          alpha_m = alpha[m]; // old value
          
          // correlation with residuals
          alpha[m] = alpha_m + arma::dot(Z.col(m), r) / Z2[m];
          
          dalpha = alpha[m] - alpha_m;
          if (std::fabs(dalpha) > 1e-20) {
            max_update = std::max(max_update, std::fabs(dalpha));
            
            // update residuals
            r -= Z.col(m) * dalpha;
          }
        }
        
        max_update = std::max(max_update, std::fabs(dalpha));
        
        // Solve lasso over ever-active set
        max_update = 0;
        L0 = 0;
        L0_het = 0;
        L1 = 0.0;
        L2 = 0.0;
        for (j = 0; j < p; j++) {
          if (in_A[j]) {
            if (grouped[j]) {
              cutoff = std::sqrt(lambda0 / (XE_scale[j] + 2*lambda2)) + 0.5*lambda1 / (XE_scale[j] + 2*lambda2);
              
              for (k=0; k<K+1; k++) {
                beta_j[k] = beta(j,k);
                beta_tilde[k] = beta(j,k);
              }
              
              for (i = 0; i < n; i++) {
                z = macc(i,j) * r[i] / n / (XE_scale[j] + 2*lambda2);
                beta_tilde[0] += z;
                
                for (k=1; k<K+1; k++) {
                  beta_tilde[k] += z * E(i,k-1);
                }
              }
              
              beta_norm = my_vecnorm(beta_tilde);
              score[j] = beta_norm;
              
              // thresholding
              if (beta_norm > cutoff) {
                // update beta
                beta.row(j) = beta_tilde * (1 - 0.5*lambda1 / beta_norm / (XE_scale[j] + 2*lambda2));
                
                L0 ++;
                L0_het ++; 
                L1 += beta_norm - 0.5*lambda1 / (XE_scale[j] + 2*lambda2);
                L2 += (beta_norm - 0.5*lambda1 / (XE_scale[j] + 2*lambda2)) * (beta_norm - 0.5*lambda1 / (XE_scale[j] + 2*lambda2));
                
                // difference in beta
                dbeta[0] = beta(j,0) - beta_j[0];
                max_dbeta = std::fabs(dbeta[0]);
                for (k=1; k<K+1; k++) {
                  dbeta[k] = beta(j,k) - beta_j[k];
                  max_dbeta = std::max(max_dbeta, std::fabs(dbeta[k]));
                }
                
              } else {
                for (k=0; k<K+1; k++) {
                  beta(j,k) = 0.0;
                  dbeta[k] = -beta_j[k];
                  max_dbeta = std::max(max_dbeta, std::fabs(beta_j[k]));
                }
                in_A[j] = false;
              }
              
            } else {
              cutoff = std::sqrt(lambda0 / (X_scale[j] + 2*lambda2)) + 0.5*lambda1 / (X_scale[j] + 2*lambda2);
              
              beta_j[0] = beta(j,0);
              beta_tilde[0] = beta(j,0);
              
              for (i = 0; i < n; i++) {
                beta_tilde[0] += macc(i,j) * r[i] / n / X_scale[j];
              }
              
              beta_norm = std::fabs(beta_tilde[0]);
              score[j] = beta_norm;
              
              // thresholding
              if (beta_norm > cutoff) {
                // update beta
                beta(j,0) = beta_tilde[0] * (1 - 0.5*lambda1 / beta_norm / (X_scale[j] + 2*lambda2));
                L0 ++;
                L1 += std::fabs(beta(j,0));
                L2 += beta(j,0) * beta(j,0);
                
                // difference in beta
                dbeta[0] = beta(j,0) - beta_j[0];
                max_dbeta = std::fabs(dbeta[0]);
                
              } else {
                beta(j,0) = 0.0;
                dbeta[0] = -beta_j[0];
                max_dbeta = std::fabs(beta_j[0]);
                in_A[j] = false;
              }
            }
            
            max_update = std::max(max_update, max_dbeta);
            
            // update residuals
            if (std::fabs(dbeta[0]) > 1e-20) {
              for (i=0; i<n; i++) {
                r(i) -= macc(i,j) * dbeta[0];
              }
            }
            if (grouped[j]) {
              for (k=1; k<K+1; k++) {
                if (std::fabs(dbeta[k]) > 1e-20) {
                  for (i=0; i<n; i++) {
                    r(i) -= macc(i,j) * E(i,k-1) * dbeta[k];
                  }
                }
              }
            }
          }
        }
        // Check for convergence
        if (max_update < thresh) {
          break;
        }
      }
      
      // Scan for violations in strong set
      violations = check_strong_set(
        in_A, in_S, macc, E, score, X_scale, XE_scale, grouped, beta, lambda0, lambda1, lambda2, r);
      if (violations == 0 || check > ncheck) {
        break;
      } else {
        check++;
      }
    }
    
    message = lambda_ix + "." + pval + "." + std::to_string(l) + ":" + std::to_string(iter[l]) + "," + std::to_string(L0) + "\n";
    Rcout << message;
    
    // objective
    metrics[l] = gLoss(r) / n;
    loss[l] = metrics[l] + lambda0*L0 + lambda1*L1 + lambda2*L2;
    nb_active[l] = L0;
    nb_interact[l] = L0_het;
    nb_checks[l] = check;
    nb_candidate[l] = n_strong;
    
    // update output
    //// SNP effects
    for (j=0; j<p; j++) {
      if (in_A[j]) {
        for (k=0; k<K+1; k++) {
          out_beta(j, (K+1)*l + k) = beta(j,k);
        }
      }
    }
    
    //// covariate effects
    for (m=0; m<M; m++) {
      out_alpha(m, l) = alpha[m];
    }
    
    // check dfmax
    if (nb_active[l] >= dfmax) {
      return List::create(_["beta"] = out_beta, 
                          _["alpha"] = out_alpha, 
                          _["loss"] = loss, 
                          _["metrics"] = metrics, 
                          _["iter"] = iter, 
                          _["messages"] = "Too many variables",
                          _["nb_active"] = nb_active, 
                          _["nb_interact"] = nb_interact, 
                          _["nb_checks"] = nb_checks, 
                          _["nb_candidate"] = nb_candidate);
    }
  }
  
  return List::create(_["beta"] = out_beta, 
                      _["alpha"] = out_alpha, 
                      _["loss"] = loss, 
                      _["metrics"] = metrics, 
                      _["iter"] = iter, 
                      _["messages"] = "Complete path",
                      _["nb_active"] = nb_active, 
                      _["nb_interact"] = nb_interact, 
                      _["nb_checks"] = nb_checks, 
                      _["nb_candidate"] = nb_candidate);
}

/******************************************************************************/

}}

#endif // #ifndef BIGSTATSR_BIGLASSO_LIN_HPP_INCLUDED
