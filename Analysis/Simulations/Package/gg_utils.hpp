#ifndef BIGSTATSR_BIGLASSO_UTILS_HPP_INCLUDED
#define BIGSTATSR_BIGLASSO_UTILS_HPP_INCLUDED

/******************************************************************************/

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using std::size_t;

/******************************************************************************/

namespace gg_bigstatsr { namespace big_gglassoUtils {


double my_vecnorm(const rowvec & beta_tilde) {
  double out=0.0;
  
  for (int k=0; k<beta_tilde.size(); k++) out += beta_tilde[k]*beta_tilde[k];
  
  out = std::sqrt(out);
  return out;
}

double my_vecnorm(const vec & beta_tilde) {
  double out=0.0;

  for (int k=0; k<beta_tilde.size(); k++) out += beta_tilde[k]*beta_tilde[k];

  out = std::sqrt(out);
  return out;
}

// separate main and interaction summaries
template <class C>
List get_summaries2(C macc, const mat& E, const vec& Y) {
    // Rcout << "gg_utils.hpp get_summaries \n";
    
    int n = macc.nrow();
    int p = macc.ncol();
    int K = E.n_cols;
    
    // Rcout << n << " " << p << " " << K << " " << M << "\n";
    
    NumericMatrix XE_scale(p, K+1); // L2 norms of main and interaction effects
    
    NumericMatrix XEY(p, K+1); // ||(XE)'Y||: correlation of main / interactions w outcome
    
    int i=0, j=0, k=0;
    
    double x=0.0, xe=0.0;
    
    for (j=0; j<p; j++) {
      if (j % 100000 == 0) {
        Rcout << j << "..\n";
      } else if (j % 10000 == 0) {
        Rcout << j << "..";
      }
      
      for (i=0; i<n; i++) {
        // SNP terms
        x = macc(i, j);
        
        // SNP-SNP
        XE_scale(j,0) += x * x / n;
        
        // SNP-Y
        XEY(j,0) += x * Y[i] / n;
        
        // interaction terms
        for (k=1; k<K+1; k++) {
          xe = x * E(i,k-1);
          
          // Interaction-Interaction
          XE_scale(j,k) += xe * xe / n;
          
          // Interaction-Y
          XEY(j,k) += xe * Y[i] / n;
          
        }
      }
    }
    
    Rcout << "\n";
    
    return List::create(_["XE_scale"] = XE_scale,
                        _["XEY"] = XEY);
}

// grouped and main-only summaries
template <class C>
List get_summaries(C macc, const mat& E, const vec& Y,
                   const int maxit = 1e3, const double tol = 1e-11) {
  // Rcout << "gg_utils.hpp get_summaries \n";
  
  int n = macc.nrow();
  int p = macc.ncol();
  int K = E.n_cols;
  
  // Rcout << n << " " << p << " " << K << " " << M << "\n";
  
  vec X_scale(p); // X'X / n
  vec XE_scale(p); // largest eigenvalue of (XE)'(XE) / n
  mat H(K+1,K+1);
  vec vec0(K+1), vec1(K+1);
  vec it(p); 
  it.fill(maxit);
  
  vec norm_XY(p); // ||(XE)'Y||: correlation of SNP w outcome
  vec XEY(K+1); // (XE)'Y: correlation of group w outcome
  vec norm_XEY(p); // ||(XE)'Y||: correlation of group w outcome
  
  double x=0.0, xe=0.0, xe0=0.0;
  double eig=0.0, eig0=0.0;
  
  // largest eigenvalue of group (utilize symmetric matrix structure of H)
  for (int j=0; j<p; j++) {
    if (j % 100000 == 0) {
      Rcout << j << "..\n";
    } else if (j % 10000 == 0) {
      Rcout << j << "..";
    }
    
    XEY.fill(0.0);
    
    H.fill(0.0);
    
    for (int i=0; i<n; i++) {
      // SNP terms
      x = macc(i, j);
      
      // SNP-SNP
      H(0,0) += x * x;
      
      // SNP-Y
      XEY[0] += x * Y[i];
      
      // interaction terms
      for (int k=1; k<K+1; k++) {
        // Interaction-Y
        XEY[k] += x * E(i,k-1) * Y[i];
        
        // SNP-Interaction
        xe = x * E(i,k-1);
        H(0,k) += x * xe;
        H(k,0) += x * xe;
        
        // Interaction-Interaction
        for (int k0=k+1; k0<K+1; k0++) {
          xe0 = x * E(i,k0-1);
          H(k,k0) += xe * xe0;
          H(k0,k) += xe * xe0;
        }
        
        // Interaction-Interaction
        H(k,k) += xe * xe;
      }
    }
    H = H / n;
    // Rcout << H << '\n';
    
    // power iterations for largest eigenvalue
    vec0.fill(1.0 / std::sqrt(K+1)); 
    eig0 = 0.0;
    eig = 0.0;
    
    for (int t=0; t<maxit; t++) {
      eig0 = eig;
      vec1 = H * vec0;
      eig = arma::dot(vec0, vec1);
      
      if (std::fabs(eig0 - eig) < tol) {
        it[j] = t;
        break;
      }
      
      // Rcout << eig0 << " " << eig << " " << std::fabs(eig0 - eig) << "\n";
      
      vec0 = vec1 / my_vecnorm(vec1);
    }
    XE_scale[j] = eig;
    X_scale[j] = H(0,0);
    
    norm_XEY[j] = my_vecnorm(XEY) / n;
    norm_XY[j] = std::fabs(XEY[0]) / n;
  }
  Rcout << "\n";
  
  return List::create(_["X_scale"] = X_scale,
                      _["XE_scale"] = XE_scale,
                      _["XY_norm"] = norm_XY,
                      _["XEY_norm"] = norm_XEY,
                      _["power_it"] = it);
}




// grouped and main-only summaries (in parallel)
template <class C>
List get_summaries_par(C macc, const mat& E, const vec& Y, const int ncores,
                   const int maxit = 1e3, const double tol = 1e-11) {
  // Rcout << "gg_utils.hpp get_summaries \n";
  
  int n = macc.nrow();
  int p = macc.ncol();
  int K = E.n_cols;
  
  vec X_scale(p); // X'X / n
  vec XE_scale(p); // largest eigenvalue of (XE)'(XE) / n
  vec it(p); 
  
  vec norm_XY(p); // ||(XE)'Y||: correlation of SNP w outcome
  vec norm_XEY(p); // ||(XE)'Y||: correlation of group w outcome
  
  // Rcout << n << " " << p << " " << K << " " << M << "\n";
  
  // largest eigenvalue of group (utilize symmetric matrix structure of H)
  int chunk_size = ceil(p / (10.0 * ncores));
  
  #pragma omp parallel num_threads(ncores)
  {
    double x=0.0, xe=0.0, xe0=0.0;
    double eig=0.0, eig0=0.0;
    
    mat H(K+1,K+1);
    vec XEY(K+1); // (XE)'Y: correlation of group w outcome
    
    vec vec0(K+1), vec1(K+1);
    it.fill(maxit);
    
    #pragma omp for schedule(dynamic, chunk_size)
    for (int j=0; j<p; j++) {
      if (j % 1000 == 0) {
      Rcout << j << "..";
      }
      
      XEY.fill(0.0);
      
      H.fill(0.0);
      
      for (int i=0; i<n; i++) {
        // SNP terms
        x = macc(i, j);
        
        // SNP-SNP
        H(0,0) += x * x;
        
        // SNP-Y
        XEY[0] += x * Y[i];
        
        // interaction terms
        for (int k=1; k<K+1; k++) {
          // Interaction-Y
          XEY[k] += x * E(i,k-1) * Y[i];
          
          // SNP-Interaction
          xe = x * E(i,k-1);
          H(0,k) += x * xe;
          H(k,0) += x * xe;
          
          // Interaction-Interaction
          for (int k0=k+1; k0<K+1; k0++) {
            xe0 = x * E(i,k0-1);
            H(k,k0) += xe * xe0;
            H(k0,k) += xe * xe0;
          }
          
          // Interaction-Interaction
          H(k,k) += xe * xe;
        }
      }
      H = H / n;
      // Rcout << H << '\n';
      
      // power iterations for largest eigenvalue
      vec0.fill(1.0 / std::sqrt(K+1)); 
      eig0 = 0.0;
      eig = 0.0;
      
      for (int t=0; t<maxit; t++) {
        eig0 = eig;
        vec1 = H * vec0;
        eig = arma::dot(vec0, vec1);
        
        if (std::fabs(eig0 - eig) < tol) {
          it[j] = t;
          break;
        }
        
        // Rcout << eig0 << " " << eig << " " << std::fabs(eig0 - eig) << "\n";
        
        vec0 = vec1 / my_vecnorm(vec1);
      }
      XE_scale[j] = eig;
      X_scale[j] = H(0,0);
      
      norm_XEY[j] = my_vecnorm(XEY) / n;
      norm_XY[j] = std::fabs(XEY[0]) / n;
    }
  }
  
  Rcout << "\n";
  
  return List::create(_["X_scale"] = X_scale,
                      _["XE_scale"] = XE_scale,
                      _["XY_norm"] = norm_XY,
                      _["XEY_norm"] = norm_XEY,
                      _["power_it"] = it);
}


/******************************************************************************/

template <class C>
vec predict(C macc_val,
            const mat& beta,
            const vec& alpha,
            const mat& Z_val,
            const mat& E_val,
            const LogicalVector& grouped,
            const LogicalVector& in_A) {

    // Rcout << "gg_utils.hpp predict \n";
    // Rcout << "rows: " << macc.nrow() << "\n";
    // Rcout << "cols: " << macc.ncol() << "\n";
    // Rcout << "Z: " << Z.n_cols << "\n";
    // Rcout << "E: " << E.n_cols << "\n";
    // Rcout << typeid(macc(0,0)).name() << "\n";
    // Rcout << macc(0,0) << "\n";
    
  int n = macc_val.nrow();
  int p = macc_val.ncol();
  int K = E_val.n_cols;
  vec pred = Z_val * alpha; // start with covariates

  // SNP effect
  for (int j = 0; j < p; j++) {
    if (in_A[j]) {
      for (int i = 0; i < n; i++) {
        // main effect
        pred[i] += macc_val(i, j) * beta(j,0);
        
        // interaction effect
        if (grouped[j]) {
          for (int k=1; k<K+1; k++) {
            pred[i] += macc_val(i, j) * E_val(i,k-1) * beta(j,k);
          }
        }
      }
    }
  }
  
  return(pred);
}
/******************************************************************************/


// check KKT conditions over features in the strong set (stronger)
template <class C>
size_t COPY_check_strong_set(LogicalVector& in_A,
                             const LogicalVector& in_S,
                             C macc,
                             const mat& E,
                             vec& score,
                             const vec& X_scale, 
                             const vec& XE_scale, 
                             const LogicalVector& grouped, 
                             const mat& beta,
                             double lambda0, double lambda1, double lambda2,
                             const vec& r) {
  
  // Rcout << "gg_utils.hpp COPY_check_strong_set \n";
  // Rcout << "rows: " << macc.nrow() << "\n";
  // Rcout << "cols: " << macc.ncol() << "\n";
  // Rcout << typeid(macc(0,0)).name() << "\n";
  // Rcout << macc(0,0) << "\n";
  
  size_t n = macc.nrow();
  size_t p = macc.ncol();
  size_t K = E.n_cols;
  size_t i, j, k, violations = 0;
  double cutoff = 0.0;
  vec beta_tilde(K+1);
  double beta_norm = 0.0;
  double z=0.0;
  
  for (j = 0; j < p; j++) {
    if (in_S[j] && !in_A[j]) {
      if (grouped[j]) {
        cutoff = std::sqrt(lambda0 / (XE_scale[j] + 2*lambda2)) + 0.5*lambda1 / (XE_scale[j] + 2*lambda2);
        
        beta_tilde.fill(0.0);
        
        for (i = 0; i < n; i++) {
          z=macc(i,j) * r[i];
          beta_tilde[0] += z;
          
          for (k=1; k<K+1; k++) {
            beta_tilde[k] += z * E(i,k-1);
          }
        }
        
        beta_norm = my_vecnorm(beta_tilde) / n / (XE_scale[j] + 2*lambda2);
      } else {
        cutoff = std::sqrt(lambda0 / (X_scale[j] + 2*lambda2)) + 0.5*lambda1 / (X_scale[j] + 2*lambda2);
        
        beta_tilde[0] = 0.0;
        
        for (i = 0; i < n; i++) {
          beta_tilde[0] += macc(i,j) * r[i];
        }
        
        beta_norm = std::fabs(beta_tilde[0]) / n / (X_scale[j] + 2*lambda2);
      }
      
      score[j] = beta_norm;
      if (score[j] > cutoff) {
        violations++;
        in_A[j] = true;
      } 
    }
  }
  
  return(violations);
}

/******************************************************************************/

template <class C>
vec predict2(C macc_val,
            const mat& beta,
            const vec& alpha,
            const mat& Z_val,
            const mat& E_val,
            const LogicalVector& grouped,
            const LogicalMatrix& in_A) {
  
  // Rcout << "gg_utils.hpp predict \n";
  // Rcout << "rows: " << macc.nrow() << "\n";
  // Rcout << "cols: " << macc.ncol() << "\n";
  // Rcout << "Z: " << Z.n_cols << "\n";
  // Rcout << "E: " << E.n_cols << "\n";
  // Rcout << typeid(macc(0,0)).name() << "\n";
  // Rcout << macc(0,0) << "\n";
  
  int n = macc_val.nrow();
  int p = macc_val.ncol();
  int K = E_val.n_cols;
  vec pred = Z_val * alpha; // start with covariates
  
  int i=0, j=0, k=0;
  
  // SNP effect
  for (j = 0; j < p; j++) {
    // main effect
    if (in_A(j,0)) {
      for (i = 0; i < n; i++) pred[i] += macc_val(i, j) * beta(j,0);
    }
    
    // interaction effect
    if (grouped[j]) {
      for (k=1; k<K+1; k++) {
        if (in_A(j,k)) {
          for (i =0; i<n; i++) pred[i] += macc_val(i, j) * E_val(i,k-1) * beta(j,k);
        }
      }
    }
  }
  
  return(pred);
}
/******************************************************************************/


// check KKT conditions over features in the strong set (stronger)
template <class C>
size_t COPY_check_strong_set2(LogicalMatrix& in_A,
                              const LogicalMatrix& in_S,
                              C macc,
                              const mat& E,
                              mat& score,
                              const mat& XE_scale, 
                              const LogicalVector& grouped, 
                              const mat& beta,
                              double lambda0, double lambda1, double lambda2,
                              const vec& r) {
  
  // Rcout << "gg_utils.hpp COPY_check_strong_set \n";
  // Rcout << "rows: " << macc.nrow() << "\n";
  // Rcout << "cols: " << macc.ncol() << "\n";
  // Rcout << typeid(macc(0,0)).name() << "\n";
  // Rcout << macc(0,0) << "\n";
  
  size_t n = macc.nrow();
  size_t p = macc.ncol();
  size_t K = E.n_cols;
  size_t i=0, j=0, k=0, violations = 0;
  double cutoff = 0.0;
  double beta_tilde=0.0;
  double z=0.0;
  
  for (j = 0; j < p; j++) {
    // Rcout << j << '.';
    // main effects
    if (in_S(j,0) && !in_A(j,0)) {
      cutoff = std::sqrt(lambda0 / (XE_scale(j,0) + 2*lambda2)) + 0.5*lambda1 / (XE_scale(j,0) + 2*lambda2);
      
      beta_tilde = 0.0;
      
      for (i = 0; i < n; i++) {
        beta_tilde += macc(i,j) * r[i];
      }
      
      score(j,0) = std::fabs(beta_tilde) / n / (XE_scale(j,0) + 2*lambda2);
      
      if (score(j,0) > cutoff) {
        violations++;
        in_A(j,0) = true;
      }
    }
    
    // interactions
    if (grouped[j]) {
      for (k=1; k<K+1; k++) {
        if (in_S(j,k) && !in_A(j,k)) {
          cutoff = std::sqrt(lambda0 / (XE_scale(j,k) + 2*lambda2)) + 0.5*lambda1 / (XE_scale(j,k) + 2*lambda2);
          
          beta_tilde = 0.0;
          
          for (i = 0; i < n; i++) {
            beta_tilde += macc(i,j) * r[i] * E(i,k-1);
          }
          
          score(j,k) = std::fabs(beta_tilde) / n / (XE_scale(j,k) + 2*lambda2);
          
          if (score(j,k) > cutoff) {
            violations++;
            in_A(j,k) = true;
          }
        }
      }
    }
  }
  
  return(violations);
}

} }

/******************************************************************************/

#endif // #ifndef BIGSTATSR_BIGLASSO_UTILS_HPP_INCLUDED
