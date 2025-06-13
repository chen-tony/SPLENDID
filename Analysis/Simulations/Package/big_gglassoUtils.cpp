// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, bigsnpr, rmio, RcppArmadillo)]]

/******************************************************************************/
#include <RcppArmadillo.h>
#include <gg_BMAcc.h>
#include <gg_BMAcc-dispatcher.h>
#include <gg_utils.hpp>

/******************************************************************************/

#define CALL_BIGSUMMARIES(ACC) {                                                \
Rcout << "my_biglassoUtils.cpp CALL_BIGSUMMARIES \n"; \
  return gg_bigstatsr::big_gglassoUtils::get_summaries(ACC, E, Y, maxit, tol);    \
}

// Dispatch function for get_summaries
// [[Rcpp::export]]
List bigsummaries(Environment BM,
                  const IntegerVector& rowInd,
                  const IntegerVector& colInd,
                  const mat& E,
                  const vec& Y,
                  const int maxit = 1e6, const double tol = 1e-11) {
  // Rcout << "big_gglassoUtils.cpp dispatch \n";
  DISPATCH_SUBMATACC(CALL_BIGSUMMARIES);
}

/******************************************************************************/

#define CALL_BIGSUMMARIES2(ACC) {                                                \
Rcout << "my_biglassoUtils.cpp CALL_BIGSUMMARIES2 \n";                           \
return gg_bigstatsr::big_gglassoUtils::get_summaries2(ACC, E, Y);    \
}

// Dispatch function for get_summaries
// [[Rcpp::export]]
List bigsummaries2(Environment BM,
                  const IntegerVector& rowInd,
                  const IntegerVector& colInd,
                  const mat& E,
                  const vec& Y) {
  // Rcout << "big_gglassoUtils.cpp dispatch \n";
  DISPATCH_SUBMATACC(CALL_BIGSUMMARIES2);
}

/******************************************************************************/

#define CALL_BIGSUMMARIES_PAR(ACC) {                                                \
Rcout << "my_biglassoUtils.cpp CALL_BIGSUMMARIES_PAR \n";                           \
return gg_bigstatsr::big_gglassoUtils::get_summaries_par(ACC, E, Y, ncores, maxit, tol);    \
}

// Dispatch function for get_summaries
// [[Rcpp::export]]
List bigsummaries_par(Environment BM,
                  const IntegerVector& rowInd,
                  const IntegerVector& colInd,
                  const mat& E,
                  const vec& Y, const int ncores,
                  const int maxit = 1e6, const double tol = 1e-11) {
  // Rcout << "big_gglassoUtils.cpp dispatch \n";
  DISPATCH_SUBMATACC(CALL_BIGSUMMARIES_PAR);
}

/******************************************************************************/
