/////////////////////////////////////////////////////////////////////////////////////
/// adapted from original bigstatsr package: https://github.com/privefl/bigstatsr ///
/////////////////////////////////////////////////////////////////////////////////////

#ifndef SLINGSHOT_UTILS_HPP_INCLUDED
#define SLINGSHOT_UTILS_HPP_INCLUDED

/******************************************************************************/

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using std::size_t;

/******************************************************************************/

namespace slingshot { namespace slingshot_utils {


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

// grouped and main-only summaries
template <class C>
List get_summaries(C macc, const mat& E, const vec& Y,
                   const int maxit = 1e3, const double tol = 1e-11) {
  
  int n = macc.nrow();
  int p = macc.ncol();
  int K = E.n_cols;
  
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
  
  int n = macc.nrow();
  int p = macc.ncol();
  int K = E.n_cols;
  
  vec X_scale(p); // X'X / n
  vec XE_scale(p); // largest eigenvalue of (XE)'(XE) / n
  vec it(p); 
  
  vec norm_XY(p); // ||(XE)'Y||: correlation of SNP w outcome
  vec norm_XEY(p); // ||(XE)'Y||: correlation of group w outcome
  
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

// make prediction
template <class C>
vec predict(C macc_tun,
            const mat& beta,
            const vec& alpha,
            const mat& Z_tun,
            const mat& E_tun,
            const LogicalVector& grouped,
            const LogicalVector& in_A) {
    
  int n = macc_tun.nrow();
  int p = macc_tun.ncol();
  int K = E_tun.n_cols;
  vec pred = Z_tun * alpha; // start with covariates

  // SNP effect
  for (int j = 0; j < p; j++) {
    if (in_A[j]) {
      for (int i = 0; i < n; i++) {
        // main effect
        pred[i] += macc_tun(i, j) * beta(j,0);
        
        // interaction effect
        if (grouped[j]) {
          for (int k=1; k<K+1; k++) {
            pred[i] += macc_tun(i, j) * E_tun(i,k-1) * beta(j,k);
          }
        }
      }
    }
  }
  
  return(pred);
}
/******************************************************************************/


// check KKT conditions over features in the strong set
template <class C>
size_t check_strong_set(LogicalVector& in_A,
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

} }

/******************************************************************************/

#endif // #ifndef BIGSTATSR_BIGLASSO_UTILS_HPP_INCLUDED
