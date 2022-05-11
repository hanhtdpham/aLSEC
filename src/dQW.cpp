#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat dQW(const arma::colvec & SS,
              const arma::colvec & RR,
              const arma::mat & U,
              const arma::mat & V,
              const arma::mat & W,
              const arma::mat & Pmk,
              const arma::rowvec & Pk,
              const IntegerMatrix & EE){
  int M = Pmk.n_rows;
  int K = W.n_rows;
  int p = W.n_cols;
  int n = U.n_rows;

  arma::mat eUiWki = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWki = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWki,0);
  arma::rowvec fvk = sum(eViWki,0);
  arma::mat gradW = -W;

  arma::mat Suk = arma::zeros(K,p);
  arma::mat Svk = arma::zeros(K,p);
  for(int i=0;i<n;i++){
    Suk = Suk + eUiWki.row(i).t()*U.row(i);
    Svk = Svk + eViWki.row(i).t()*V.row(i);
  }

  for(int k=0; k<K; k++){

    for(int m=0; m<M; m++){
      gradW.row(k) = gradW.row(k) +
        Pmk(m,k)*(
            U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1) -
            Suk.row(k)/fuk(k) -
            (Svk.row(k) - eViWki(EE(m,0) - 1,k)*V.row(EE(m,0) - 1))/
              (fvk(k) - eViWki(EE(m,0) - 1,k))
        );
    }

  }

  return gradW;
}
