#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat dQUV(const arma::colvec & SS,
               const arma::colvec & RR,
               const arma::mat & U,
               const arma::mat & V,
               const arma::mat & W,
               const arma::mat & Pmk,
               const arma::rowvec & Pk,
               const arma::mat & Pmki1,
               const arma::mat & Pmki2,
               const double & tauS,
               const double & tauR,
               const double & tauU,
               const double & tauV,
               const IntegerMatrix & EE,
               const IntegerVector & Mi1Index,
               const IntegerMatrix & Mi1){
  int M = EE.nrow();
  int K = W.n_rows;
  int p = U.n_cols;
  int n = U.n_rows;

  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  arma::colvec hk = arma::zeros(K);

  for(int m=0; m<M; m++){
    for(int k=0; k<K; k++){
      hk(k) = hk(k) +
        Pmk(m,k)/(fvk(k) - eUiWk(EE(m,0) - 1,k));
    }
  }

  arma::mat gradUV = arma::zeros(n,2*(p + 1));
  arma::mat WTilde = arma::ones(K,p + 1);
  WTilde.cols(1,p) = W;

  for(int i=0;i<n;i++){

    for(int k=0; k<K; k++){
      gradUV.row(i).cols(0,p) = gradUV.row(i).cols(0,p) +
        (Pmki1(i,k) - Pk(k)/fuk(k)*eUiWk(i,k))*WTilde.row(k);
      gradUV.row(i).cols(p + 1,2*(p + 1) - 1) = gradUV.row(i).cols(p + 1,2*(p + 1) - 1) +
        (Pmki2(i,k) - eViWk(i,k)*(hk(k) - Pmki1(i,k)/(fvk(k) - eViWk(i,k))) )*WTilde.row(k);
    }

    gradUV(i,0) = gradUV(i,0) - SS(i)*tauS;
    gradUV.row(i).cols(1,p) = gradUV.row(i).cols(1,p) - U.row(i)*tauU;
    gradUV(i,p + 1) = gradUV(i,p + 1) - RR(i)*tauR;
    gradUV.row(i).cols(p + 2,2*(p + 1) - 1) = gradUV.row(i).cols(p + 2,2*(p + 1) - 1) - V.row(i)*tauV;

  }

  return gradUV;
}
