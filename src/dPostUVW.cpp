#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List dPostUVW(const arma::colvec & SS,
                    const arma::colvec & RR,
                    const arma::mat & U,
                    const arma::mat & V,
                    const arma::mat & W,
                    const IntegerVector & z,
                    const double & tauS,
                    const double & tauR,
                    const double & tauU,
                    const double & tauV,
                    const IntegerMatrix & EE,
                    const arma::colvec & n_k,
                    const IntegerVector & Mi1Index,
                    const IntegerMatrix & Mi1,
                    const IntegerVector & Mi2Index,
                    const IntegerMatrix & Mi2){
  int M = EE.nrow();
  int K = W.n_rows;
  int n = U.n_rows;
  int p = U.n_cols;
  
  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  arma::mat gradUV = arma::zeros(n,2*(p + 1));
  arma::mat gradW = arma::zeros(K,p);
  
  double helper = 0;
  arma::cube Suk_elements = arma::zeros(K,p,n);
  arma::mat Suk = arma::zeros(K,p);
  arma::cube Svk_elements = arma::zeros(K,p,n);
  arma::mat Svk = arma::zeros(K,p);
  arma::colvec s_tilde = arma::zeros(K);
  
  
  // for(int k=0;k<K;k++){
  //   for(int m=0;m<M;m++){
  //     s_tilde(k) = s_tilde(k) +
  //       1/(fvk(k) - eViWk(EE(m,0) - 1,k));
  //   }
  // }
  for(int m=0;m<M;m++){
    s_tilde(z(m) - 1) = s_tilde(z(m) - 1) +
      1/(fvk(z(m) - 1) - eViWk(EE(m,0) - 1,z(m) - 1));
  }
  
  for(int i=0;i<n;i++){
    
    // (S,U)
    if(Mi1Index(i)>0){
      if(Mi1Index(i) == 1){
        gradUV(i,0) = 1.0;
        gradUV.submat(i,1,i,p) = W.row(z(Mi1(i,0) - 1) - 1);
      }else{
        gradUV(i,0) = 1.0*Mi1Index(i);
        for(int j=0;j<Mi1Index(i);j++){
          gradUV.submat(i,1,i,p) = gradUV.submat(i,1,i,p) +
            W.row(z(Mi1(i,j) - 1) - 1);
        }
      }
    }
    for(int k=0; k<K; k++){
      helper = n_k(k)*eUiWk(i,k)/fuk(k);
      gradUV(i,0) = gradUV(i,0) - helper;
      gradUV.submat(i,1,i,p) = gradUV.submat(i,1,i,p) -
        helper*W.row(k);
    }
    gradUV(i,0) = gradUV(i,0) - tauS*SS(i);
    gradUV.submat(i,1,i,p) = gradUV.submat(i,1,i,p) -
      tauU*U.row(i);
    
    // (R,V)
    if(Mi2Index(i)>0){
      if(Mi2Index(i) == 1){
        gradUV(i,p + 1) = 1.0;
        gradUV.submat(i,p + 2,i,2*p + 1) = W.row(z(Mi2(i,0) - 1) - 1);
      }else{
        gradUV(i,p + 1) = 1.0*Mi2Index(i);
        for(int j=0;j<Mi2Index(i);j++){
          gradUV.submat(i,p + 2,i,2*p + 1) = gradUV.submat(i,p + 2,i,2*p + 1) +
            W.row(z(Mi2(i,j) - 1) - 1);
        }
      }
    } 
    
    for(int k=0; k<K; k++){
      helper = 0;
      if(Mi1Index(i)>0){
        for(int j=0;j<Mi1Index(i);j++){
          if(z(Mi1(i,j) - 1) - 1 == k){
            helper += 1.0;
          }
        }
      }
      helper = eViWk(i,k)*(s_tilde(k) - helper/(fvk(k) - eViWk(i,k)));
      gradUV(i,p + 1) = gradUV(i,p + 1) - helper;
      gradUV.submat(i,p + 2,i,2*p + 1) = gradUV.submat(i,p + 2,i,2*p + 1) -
        helper*W.row(k);
    }
    gradUV(i,p + 1) = gradUV(i,p + 1) - tauR*RR(i);
    gradUV.submat(i,p + 2,i,2*p + 1) = gradUV.submat(i,p + 2,i,2*p + 1) -
      tauV*V.row(i);
  }
  
  // W
  for(int i=0;i<n;i++){
    Suk_elements.slice(i) = eUiWk.row(i).t()*U.row(i);
    Suk = Suk + Suk_elements.slice(i);
    Svk_elements.slice(i) = eViWk.row(i).t()*V.row(i);
    Svk = Svk + Svk_elements.slice(i);
  }
  for(int m=0;m<M;m++){
    gradW.row(z(m) - 1) = gradW.row(z(m) - 1) +
      U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1) - 
      Suk.row(z(m) - 1)/fuk(z(m) - 1) -
      (Svk.row(z(m) - 1) - Svk_elements.slice(EE(m,0) - 1).row(z(m) - 1))/
        (fvk(z(m) - 1) - eViWk(EE(m,0) - 1,z(m) - 1));
  }
  gradW = gradW - W;
  
  return List::create(Named("UV") = gradUV,
                      Named("W") = gradW); //,Named("test") = helper);
}
