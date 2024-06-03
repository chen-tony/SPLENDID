// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, bigsnpr, rmio, RcppArmadillo)]]

/////////////////////////////////////////////////////////////////////////////////////
/// adapted from original bigstatsr package: https://github.com/privefl/bigstatsr ///
/////////////////////////////////////////////////////////////////////////////////////

/******************************************************************************/
#include <RcppArmadillo.h>
#include <slingshot_BMAcc.h>
#include <slingshot_BMAcc-dispatcher.h>
#include <utils.hpp>

/******************************************************************************/

#define CALL_BIGSUMMARIES(ACC) {                                                \
  return slingshot::slingshot_utils::get_summaries(ACC, E, Y, maxit, tol);    \
}

// Dispatch function for get_summaries
// [[Rcpp::export]]
List bigsummaries(Environment BM,
                  const IntegerVector& rowInd,
                  const IntegerVector& colInd,
                  const mat& E,
                  const vec& Y,
                  const int maxit = 1e6, const double tol = 1e-11) {
  DISPATCH_SUBMATACC(CALL_BIGSUMMARIES);
}

/******************************************************************************/

#define CALL_BIGSUMMARIES_PAR(ACC) {                                                \
return slingshot::slingshot_utils::get_summaries_par(ACC, E, Y, ncores, maxit, tol);    \
}

// Dispatch function for get_summaries
// [[Rcpp::export]]
List bigsummaries_par(Environment BM,
                  const IntegerVector& rowInd,
                  const IntegerVector& colInd,
                  const mat& E,
                  const vec& Y, const int ncores,
                  const int maxit = 1e6, const double tol = 1e-11) {
  DISPATCH_SUBMATACC(CALL_BIGSUMMARIES_PAR);
}

/******************************************************************************/
