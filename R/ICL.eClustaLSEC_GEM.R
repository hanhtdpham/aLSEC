#' ICL for edge clustering model using a finite mixture prior
#'
#' Compute the Integrated Complete Likelihood for edge clustering object
#'
#' Higher is better.
#'
#' @param eCl object of class 'eClustaLSEC_GEM'
#' @param A igraph object
#' @param ... Further parameters to use in \code{optim()}.
ICL.eClustaLSEC_GEM = function(eCl,A,...){
  n = vcount(A)
  K = nrow(eCl$estimates$W)
  p = ncol(eCl$estimates$U)
  M = ecount(A)

  zHat = apply(eCl$estimates$Pmk,1,which.max)
  n_k = sapply(1:K,function(k) sum(zHat == k))

  EE = cbind(ends(A,1:ecount(A),FALSE))
  temp = indexEdges(EE,n)
  Mi1 = temp$Mi1[,1:max(temp$Mi1Index)]
  Mi2 = temp$Mi2[,1:max(temp$Mi2Index)]
  Mi1Index = temp$Mi1Index
  Mi2Index = temp$Mi2Index
  rm(temp)

  wrappers = list()
  wrappers$vec_grad = function(obj){
    c(obj$UV[,1],obj$UV[,p+2],c(obj$UV[,1 + 1:p]),c(obj$UV[,p + 2 + 1:p]),c(obj$W))
  }
  wrappers$negCondLik = function(x){
    SS = x[1:n]
    RR = x[n + 1:n]
    UU = matrix(x[2*n + 1:(n*p)],n,p)
    VV = matrix(x[2*n + n*p + 1:(n*p)],n,p)
    WW = matrix(x[2*n + 2*n*p + 1:(K*p)],K,p)

    return(-evalConditionalLik(zHat,SS,RR,UU,VV,WW,EE))
  }
  wrappers$grad_negCondLik = function(x){
    SS = x[1:n]
    RR = x[n + 1:n]
    UU = matrix(x[2*n + 1:(n*p)],n,p)
    VV = matrix(x[2*n + n*p + 1:(n*p)],n,p)
    WW = matrix(x[2*n + 2*n*p + 1:(K*p)],K,p)

    grad_0 =
      dCondLik(SS,RR,UU,VV,WW,zHat,EE,n_k,
               Mi1Index,Mi1,Mi2Index,Mi2)
    return(-wrappers$vec_grad(grad_0))
  }

  Opt = optim(par=c(eCl$estimates$S,eCl$estimates$R,
                    c(eCl$estimates$U),c(eCl$estimates$V),
                    c(eCl$estimates$W)),
              fn = wrappers$negCondLik,
              gr = wrappers$grad_negCondLik,
              method="CG",...)

  ICL = -Opt$value -
    log(M)*( 2^(is_directed(A))*(vcount(A)+length(eCl$estimates$U)) + length(eCl$estimates$W)) +
    lgamma(K*eCl$estimates$alpha) + sum(lgamma(n_k +eCl$estimates$alpha)) -
    K*lgamma(eCl$estimates$alpha) - lgamma(M + K*eCl$estimates$alpha)

  return(c(ICL = ICL))
}
