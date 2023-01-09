#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double evalConditionalLik(const IntegerVector & z,
                          const arma::colvec & n_k,
                          double & lambda,
                          const arma::colvec & SS,
                          const arma::colvec & RR,
                          const arma::mat & U,
                          const arma::mat & V,
                          const arma::mat & W,
                          const IntegerMatrix & EE){
  int M = EE.nrow();
  int K = W.n_rows;
  // int n = U.n_rows;
  // int p = U.n_cols;

  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);

  double ret = 0;

  for(int m = 0;m<M;m++){
    ret +=
      SS(EE(m,0) - 1) + RR(EE(m,1) - 1) + as_scalar(
          (U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1))*W.row(z(m) - 1).t() -
            log(fuk(z(m) - 1)) - log(fvk(z(m) - 1) - eViWk(EE(m,0) - 1,z(m) - 1)) );
  }

  ret +=
    lgamma(K*exp(lambda)) + sum(lgamma(n_k+exp(lambda))) -
    K*lgamma(exp(lambda)) - lgamma(M + K*exp(lambda)) + lambda;

  return ret;
}
