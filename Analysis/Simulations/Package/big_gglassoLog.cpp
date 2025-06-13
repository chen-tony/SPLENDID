// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, bigsnpr, rmio, RcppArmadillo)]]

/******************************************************************************/
#include <gg_BMAcc.h>
#include <gg_BMAcc-dispatcher.h>
#include <gg_logistic.hpp>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

/******************************************************************************/

#define CALL_COPY_CDFIT_BINOMIAL_HSR(ACC, ACC_val) {                                                 \
Rcout << "gg_biglassoLog.cpp CALL_COPY_CDFIT_BINOMIAL_HSR \n";                                       \
return gg_bigstatsr::big_gglassoLog::COPY_cdfit_binomial_hsr(ACC, Y, Z, E, Y0, alpha0,                          \
                                                             ACC_val, Y_val, Z_val, E_val,           \
                                                             lambda, Z2, X_scale, XE_scale, XY_norm, XEY_norm, grouped, \
                                                             tol, maxit, dfmax,                      \
                                                             n_abort, nlam_min, nsupp, ncheck);             \
}

// Dispatch function for COPY_cdfit_gaussian_hsr
// [[Rcpp::export]]
List COPY_cdfit_binomial_hsr(Environment BM,
                             const vec& Y,
                             const IntegerVector& rowInd,
                             const IntegerVector& colInd,
                             const mat& Z,
                             const mat& E,
                             const vec& Y0,
                             const vec& alpha0,
                             const vec& Y_val,
                             const IntegerVector& rowInd_val,
                             const mat& Z_val,
                             const mat& E_val,
                             const mat& lambda,
                             const vec& Z2,
                             const vec& X_scale,
                             const vec& XE_scale,
                             const vec& XY_norm,
                             const vec& XEY_norm,
                             const LogicalVector& grouped,
                             double tol,
                             int maxit,
                             int dfmax,
                             int n_abort,
                             int nlam_min,
                             int nsupp,
                             int ncheck) {
  
  Rcout << "big_gglassoLog.cpp dispatch \n";
  
  DISPATCH_SUBMATACC_VAL(CALL_COPY_CDFIT_BINOMIAL_HSR)
}

/******************************************************************************/
