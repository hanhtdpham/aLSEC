#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector drawZ(const arma::colvec & SS,
                          const arma::colvec & RR,
                          const arma::mat & U,
                          const arma::mat & V,
                          const arma::mat & W,
                          const arma::rowvec & pi,
                          const IntegerMatrix & EE){
  int M = EE.nrow();
  int K = W.n_rows;

  IntegerVector z(M);
  IntegerVector KVec(K);
  for(int k=0;k<K;k++){
    KVec(k) = k + 1;
  }

  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);

  Rcpp::NumericVector Pmk(K);

  for(int m = 0;m<M;m++){
    Pmk =
      eUiWk.row(EE(m,0) - 1)%eViWk.row(EE(m,1) - 1)%pi/fuk/(fvk - eViWk.row(EE(m,0) - 1));
    Pmk = Pmk/sum(Pmk);
    z(m) = Rcpp::as<double>(RcppArmadillo::sample(KVec, 1, FALSE, Pmk));
  }

  return z;
}
