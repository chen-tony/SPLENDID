#ifndef BM_ACC_DISPATCH_H
#define BM_ACC_DISPATCH_H

/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

#define DISPATCH_ACC(CALL, ACC, RAWACC) {                                      \
Rcout << "my DISPATCH_ACC \n"; \
                                                                               \
  XPtr<FBM> xpBM = BM["address"];                                              \
                                                                               \
  if (BM.exists("code256")) {                                                  \
Rcout << 256 << '\n'; \
    CALL(RAWACC);                                                              \
  } else {                                                                     \
    switch(xpBM->matrix_type()) {                                              \
    case 8:                                                                    \
Rcout << 8 << '\n'; \
      CALL(ACC(double))                                                        \
    case 4:                                                                    \
Rcout << 4 << '\n'; \
      CALL(ACC(int))                                                           \
    case 6:                                                                    \
Rcout << 6 << '\n'; \
      CALL(ACC(float))                                                         \
    case 1:                                                                    \
Rcout << 1 << '\n'; \
      CALL(ACC(unsigned char))                                                 \
    case 2:                                                                    \
Rcout << 2 << '\n'; \
      CALL(ACC(unsigned short))                                                \
    default:                                                                   \
      throw Rcpp::exception(ERROR_TYPE);                                       \
    }                                                                          \
  }                                                                            \
}

/******************************************************************************/

#define MATACC(T) BMAcc<T>(xpBM)
#define RAWMATACC BMCode256Acc(xpBM, BM["code256"])

#define DISPATCH_MATACC(CALL) { \
Rcout << "my DISPATCH_MATACC \n"; \
DISPATCH_ACC(CALL, MATACC, RAWMATACC); \
}

#define SUBMATACC(T) SubBMAcc<T>(xpBM, rowInd, colInd, 1)
#define SUBMATACC_VAL(T) SubBMAcc<T>(xpBM, rowInd_val, colInd, 1)

#define RAWSUBMATACC SubBMCode256Acc(xpBM, rowInd, colInd, BM["code256"], 1)
#define RAWSUBMATACC_VAL SubBMCode256Acc(xpBM, rowInd_val, colInd, BM["code256"], 1)

#define DISPATCH_SUBMATACC(CALL) { \
Rcout << "my DISPATCH_SUBMATACC \n"; \
DISPATCH_ACC(CALL, SUBMATACC, RAWSUBMATACC); \
}

#define DISPATCH_SUBMATACC_VAL(CALL) {                                      \
Rcout << "my DISPATCH_SUBMATACC_VAL \n"; \
                                                                               \
  XPtr<FBM> xpBM = BM["address"];                                              \
                                                                               \
  if (BM.exists("code256")) {                                                  \
Rcout << 256 << '\n'; \
    CALL(RAWSUBMATACC, RAWSUBMATACC_VAL);                                \
  } else {                                                                     \
    switch(xpBM->matrix_type()) {                                              \
    case 8:                                                                    \
Rcout << 8 << '\n'; \
      CALL(SUBMATACC(double), SUBMATACC_VAL(double))             \
    case 4:                                                                    \
Rcout << 4 << '\n'; \
      CALL(SUBMATACC(int), SUBMATACC_VAL(int))                \
    case 6:                                                                    \
Rcout << 6 << '\n'; \
      CALL(SUBMATACC(float), SUBMATACC_VAL(float))              \
    case 1:                                                                    \
Rcout << 1 << '\n'; \
      CALL(SUBMATACC(unsigned char),  SUBMATACC_VAL(unsigned char))      \
    case 2:                                                                    \
Rcout << 2 << '\n'; \
      CALL(SUBMATACC(unsigned short), SUBMATACC_VAL(unsigned short))     \
    default:                                                                   \
      throw Rcpp::exception(ERROR_TYPE);                                       \
    }                                                                          \
  }                                                                            \
}

/******************************************************************************/

#define SUBMATCOVACC(T) SubMatCovAcc<T>(xpBM, row_idx,     col_idx, covar,     1)

#define SUBMATCOVACC_VAL(T) SubMatCovAcc<T>(xpBM, row_idx_val, col_idx, covar_val, 1)

#define RAWSUBMATCOVACC     RawSubMatCovAcc(xpBM, row_idx,     col_idx, covar,     BM["code256"], 1)

#define RAWSUBMATCOVACC_VAL RawSubMatCovAcc(xpBM, row_idx_val, col_idx, covar_val, BM["code256"], 1)

#define DISPATCH_SUBMATCOVACC(CALL)  {                                      \
Rcout << "DISPATCH_SUBMATCOVACC \n"; \
DISPATCH_ACC(CALL, SUBMATCOVACC, RAWSUBMATCOVACC);\
}

#define DISPATCH_SUBMATCOVACC_VAL(CALL) {                                      \
Rcout << "DISPATCH_SUBMATCOVACC_VAL \n"; \
                                                                               \
  XPtr<FBM> xpBM = BM["address"];                                              \
                                                                               \
  if (BM.exists("code256")) {                                                  \
Rcout << 256 << '\n'; \
    CALL(RAWSUBMATCOVACC, RAWSUBMATCOVACC_VAL);                                \
  } else {                                                                     \
    switch(xpBM->matrix_type()) {                                              \
    case 8:                                                                    \
Rcout << 8 << '\n'; \
      CALL(SUBMATCOVACC(double),         SUBMATCOVACC_VAL(double))             \
    case 4:                                                                    \
Rcout << 4 << '\n'; \
      CALL(SUBMATCOVACC(int),            SUBMATCOVACC_VAL(int))                \
    case 6:                                                                    \
Rcout << 6 << '\n'; \
      CALL(SUBMATCOVACC(float),          SUBMATCOVACC_VAL(float))              \
    case 1:                                                                    \
Rcout << 1 << '\n'; \
      CALL(SUBMATCOVACC(unsigned char),  SUBMATCOVACC_VAL(unsigned char))      \
    case 2:                                                                    \
Rcout << 2 << '\n'; \
      CALL(SUBMATCOVACC(unsigned short), SUBMATCOVACC_VAL(unsigned short))     \
    default:                                                                   \
      throw Rcpp::exception(ERROR_TYPE);                                       \
    }                                                                          \
  }                                                                            \
}

/******************************************************************************/

#endif // BM_ACC_DISPATCH_H
