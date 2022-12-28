#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double evalConditionalPost(const IntegerVector & z,
                           const arma::colvec & SS,
                           const arma::colvec & RR,
                           const arma::mat & U,
                           const arma::mat & V,
                           const arma::mat & W,
                           const double & tauS,
                           const double & tauR,
                           const double & tauU,
                           const double & tauV,
                           const IntegerMatrix & EE,
                           const arma::colvec & n_k,
                           const double & lambda,
                           const double & a_a,
                           const double & b_a){
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

  ret = ret - 0.5*tauS*sum(SS%SS) - 0.5*tauR*sum(RR%RR) -
    0.5*tauU*accu(U%U) - 0.5*tauV*accu(V%V) - 0.5*accu(W%W) +
    lgamma(K*exp(lambda)) - K*lgamma(exp(lambda)) +
    sum(lgamma(n_k + exp(lambda))) - lgamma(K*exp(lambda) + M) +
    a_a*lambda - b_a*exp(lambda);

  return ret;
}
