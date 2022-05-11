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
double QW(const arma::colvec & SS,
          const arma::colvec & RR,
          const arma::mat & U,
          const arma::mat & V,
          const arma::mat & W,
          const arma::mat & Pmk,
          const arma::rowvec & Pk,
          const IntegerMatrix & EE){
  int M = Pmk.n_rows;
  int K = W.n_rows;

  arma::mat eUiWki = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWki = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWki,0);
  arma::rowvec fvk = sum(eViWki,0);
  double ret = -0.5*accu(W%W);

  # pragma omp parallel for collapse(2)
  for(int k=0; k<K; k++){
    for(int m=0; m<M; m++){
      # pragma omp atomic
      ret += Pmk(m,k)*(
        dot(U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1), W.row(k)) -
        log(fuk(k)) -
        log(fvk(k) - eViWki(EE(m,0) - 1,k))
      );
    }
  }


  return ret;
}
