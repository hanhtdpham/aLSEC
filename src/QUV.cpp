#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#if defined(_OPENMP)
#include <omp.h>
#else
#pragma message("Warning: OpenMP is not available")
// Define functions we use which are usually found in omp.h
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_num_procs() { return 1; }
inline void omp_set_num_threads(int nthread) {}
inline void omp_set_dynamic(int flag) {}
#endif

using namespace Rcpp;

// [[Rcpp::export]]
double QUV(const arma::colvec & SS,
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
           const IntegerMatrix & EE){
  int M = EE.nrow();
  int K = W.n_rows;

  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);

  double ret = 0;

  # pragma omp parallel for collapse(2)
  for(int k=0; k<K; k++){
    for(int m = 0;m<M;m++){
      # pragma omp atomic
      ret += Pmk(m,k)*( SS(EE(m,0) - 1) + RR(EE(m,1) - 1) +
        dot(U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1), W.row(k)) -
        log(fuk(k)) - log( fvk(k) - eViWk(EE(m,0) - 1,k) ) );
    }
  }

  ret += -
    0.5*accu(U%U)*tauU - 0.5*accu(V%V)*tauV -
    0.5*sum(SS%SS)*tauS - 0.5*sum(RR%RR)*tauR;

  return ret;
}
