// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, bigsnpr, rmio, RcppArmadillo)]]

/////////////////////////////////////////////////////////////////////////////////////
/// adapted from original bigstatsr package: https://github.com/privefl/bigstatsr ///
/////////////////////////////////////////////////////////////////////////////////////

/******************************************************************************/
#include <splendid_BMAcc.h>
#include <splendid_BMAcc-dispatcher.h>
#include <linear.hpp>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

/******************************************************************************/


#define CALL_CDFIT_GAUSSIAN_PAR(ACC, ACC_tun) {                                                                                                                                                                 \
return splendid::splendid_linear::cdfit_gaussian_par(ACC, y, Z, E,                                                                \
                                                                 ACC_tun, y_tun, Z_tun, E_tun,                                                \
                                                                 lambda, Z2, X_scale, XE_scale, XY_norm, XEY_norm, grouped,                   \
                                                                 lambda_ix, pval, tol, maxit, dfmax,                                          \
                                                                 n_abort, nlam_min, ncheck);                                                  \
}

// Dispatch function for cdfit_gaussian_par
// [[Rcpp::export]]
List cdfit_gaussian_par(Environment BM,
                                 const vec& y,
                                 const IntegerVector& rowInd,
                                 const IntegerVector& colInd,
                                 const mat& Z,
                                 const mat& E,
                                 const vec& y_tun,
                                 const IntegerVector& rowInd_val,
                                 const mat& Z_tun,
                                 const mat& E_tun,
                                 const mat& lambda,
                                 const vec& Z2,
                                 const vec& X_scale,
                                 const vec& XE_scale,
                                 const vec& XY_norm,
                                 const vec& XEY_norm,
                                 const LogicalVector& grouped,
                                 const std::string & lambda_ix,
                                 const std::string & pval,
                                 double tol,
                                 int maxit,
                                 int dfmax,
                                 int n_abort,
                                 int nlam_min,
                                 int ncheck) {
  
  DISPATCH_SUBMATACC_VAL(CALL_CDFIT_GAUSSIAN_PAR)
}

/******************************************************************************/

#define CALL_CDFIT_GAUSSIAN_PAR0(ACC) {                                                                                                                                                                          \
return splendid::splendid_linear::cdfit_gaussian_par0(ACC, y, Z, E,                                                                \
                                                          lambda, Z2, X_scale, XE_scale, XY_norm, XEY_norm, grouped,                   \
                                                          lambda_ix, pval, tol, maxit, dfmax,                                          \
                                                          ncheck);                                                                     \
}

// Dispatch function for cdfit_gaussian_par0
// [[Rcpp::export]]
List cdfit_gaussian_par0(Environment BM,
                                  const vec& y,
                                  const IntegerVector& rowInd,
                                  const IntegerVector& colInd,
                                  const mat& Z,
                                  const mat& E,
                                  const mat& lambda,
                                  const vec& Z2,
                                  const vec& X_scale,
                                  const vec& XE_scale,
                                  const vec& XY_norm,
                                  const vec& XEY_norm,
                                  const LogicalVector& grouped,
                                  const std::string & lambda_ix,
                                  const std::string & pval,
                                  double tol,
                                  int maxit,
                                  int dfmax,
                                  int ncheck) {
  
  DISPATCH_SUBMATACC(CALL_CDFIT_GAUSSIAN_PAR0)
}

/******************************************************************************/
