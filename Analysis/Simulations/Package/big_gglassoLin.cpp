// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, bigsnpr, rmio, RcppArmadillo)]]

/******************************************************************************/
#include <gg_BMAcc.h>
#include <gg_BMAcc-dispatcher.h>
#include <gg_linear.hpp>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

/******************************************************************************/


#define CALL_COPY_CDFIT_GAUSSIAN_HSR(ACC, ACC_val) {                                                                    \
Rcout << "gg_biglassoLin.cpp CALL_COPY_CDFIT_GAUSSIAN_HSR \n";                                                           \
return gg_bigstatsr::big_gglassoLin::COPY_cdfit_gaussian_hsr(ACC, y, Z, E,                                              \
                                                              ACC_val, y_val, Z_val, E_val,                              \
                                                              lambda, Z2, X_scale, XE_scale, XY_norm, XEY_norm, grouped, \
                                                              tol, maxit, dfmax,                                         \
                                                              n_abort, nlam_min, ncheck);                                \
}

// Dispatch function for COPY_cdfit_gaussian_hsr
// [[Rcpp::export]]
List COPY_cdfit_gaussian_hsr(Environment BM,
                              const vec& y,
                              const IntegerVector& rowInd,
                              const IntegerVector& colInd,
                              const mat& Z,
                              const mat& E,
                              const vec& y_val,
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
                              int ncheck) {
  
  Rcout << "big_gglassoLin.cpp dispatch \n"; 
  
  DISPATCH_SUBMATACC_VAL(CALL_COPY_CDFIT_GAUSSIAN_HSR)
}

/******************************************************************************/

#define CALL_COPY_CDFIT_GAUSSIAN_HSR0(ACC) {                                                                     \
Rcout << "gg_biglassoLin.cpp CALL_COPY_CDFIT_GAUSSIAN_HSR0 \n";                                                           \
return gg_bigstatsr::big_gglassoLin::COPY_cdfit_gaussian_hsr0(ACC, y, Z, E,                                               \
                                                             lambda, Z2, X_scale, XE_scale, XY_norm, XEY_norm, grouped,  \
                                                             tol, maxit, dfmax,                                          \
                                                             ncheck);                                 \
}

// Dispatch function for COPY_cdfit_gaussian_hsr
// [[Rcpp::export]]
List COPY_cdfit_gaussian_hsr0(Environment BM,
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
                              double tol,
                              int maxit,
                              int dfmax,
                              int ncheck) {
  
  Rcout << "big_gglassoLin.cpp dispatch \n";
  
  DISPATCH_SUBMATACC(CALL_COPY_CDFIT_GAUSSIAN_HSR0)
}

/******************************************************************************/

#define CALL_COPY_CDFIT_GAUSSIAN_HSR0_2(ACC) {                                                                              \
Rcout << "gg_biglassoLin.cpp CALL_COPY_CDFIT_GAUSSIAN_HSR0 \n";                                                           \
return gg_bigstatsr::big_gglassoLin::COPY_cdfit_gaussian_hsr0_2(ACC, y, Z, E,                                               \
                                                              lambda, Z2, XE_scale, XEY, grouped,  \
                                                              tol, maxit, dfmax,                                          \
                                                              ncheck);                                                    \
}

// Dispatch function for COPY_cdfit_gaussian_hsr
// [[Rcpp::export]]
List COPY_cdfit_gaussian_hsr0_2(Environment BM,
                              const vec& y,
                              const IntegerVector& rowInd,
                              const IntegerVector& colInd,
                              const mat& Z,
                              const mat& E,
                              const mat& lambda,
                              const vec& Z2,
                              const mat& XE_scale,
                              const mat& XEY,
                              const LogicalVector& grouped,
                              double tol,
                              int maxit,
                              int dfmax,
                              int ncheck) {
  
  Rcout << "big_gglassoLin.cpp dispatch \n";
  
  DISPATCH_SUBMATACC(CALL_COPY_CDFIT_GAUSSIAN_HSR0_2)
}

/******************************************************************************/


