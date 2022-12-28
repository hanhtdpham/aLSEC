#' Summarize HMC samples
#'
#' 'summarize_eClustaLSEC' summarizes the HMC output for the LSEC model.
#' NOTE: \code{eClustaLSEC_PostProcess} should be used first!
#'
#' @param HMC_obj Object of class 'eClustaLSEC_HMC'.
#' @param burnin integer. Number of leading posterior draws to remove.
#' @param thin integer. Only each \code{Thin}-th sample will be processed.
#' @param summary_function function, e.g. mean or median. Function used to summarize the
#' posterior samples.
#' @return Object of class 'eClustaLSEC_HMC_summary'
summarize_eClustaLSEC = function(HMC_obj,
                                 burnin=1,
                                 thin=1,
                                 summary_function=mean){
  nSims = dim(HMC_obj$U)[1]

  iter_seq = seq(burnin,nSims,by=thin)

  HMC_summary = list()
  HMC_summary$z = matrix(0,ncol(HMC_obj$z),HMC_obj$K)
  for(i in 1:ncol(HMC_obj$z)){
    for(k in 1:HMC_obj$K){
      HMC_summary$z[i,k] = mean(HMC_obj$z[,i] == k)
    }
  }
  HMC_summary$S = apply(HMC_obj$S[iter_seq,],2,summary_function)
  HMC_summary$R = apply(HMC_obj$R[iter_seq,],2,summary_function)
  HMC_summary$U = apply(HMC_obj$U[iter_seq,,],2:3,summary_function)
  HMC_summary$V = apply(HMC_obj$V[iter_seq,,],2:3,summary_function)
  HMC_summary$W = apply(HMC_obj$W[iter_seq,,],2:3,summary_function)
  HMC_summary$pi = apply(HMC_obj$pi[iter_seq,],2,summary_function)
  HMC_summary$tauS = summary_function(HMC_obj$tauS[iter_seq])
  HMC_summary$tauR = summary_function(HMC_obj$tauR[iter_seq])
  HMC_summary$tauU = summary_function(HMC_obj$tauU[iter_seq])
  HMC_summary$tauV = summary_function(HMC_obj$tauV[iter_seq])
  HMC_summary$alpha = summary_function(HMC_obj$alpha[iter_seq])

  class(HMC_summary) = "eClustaLSEC_HMC_summary"
  return(HMC_summary)
}
