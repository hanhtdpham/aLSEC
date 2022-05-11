#' Automated edge clustering using a finite mixture prior
#'
#' 'eClustaLSEC_GEM' implements the variational Bayes generalized EM algorithm for clustering edges of a
#' network using a latent space approach. Use multiple starting points to find the best mode!
#'
#' @param A igraph object
#' @param K integer. Maximum number of clusters
#' @param p integer. Number of latent dimensions
#' @param maxIter integer. Maximum number of iterations of the GEM algorithm
#' @param maxIterVB integer. Maximum number of iterations of the VB algorithm for approximating E step
#' @param eps numeric. Convergence threshold for relative change of ELBO in variational E step
#' @param CGSteps integer.  Maximum number of conjugate gradient steps during each
#' iteration for U and V
#' @param QNSteps integer.  Maximum number of quasi-Newton steps during each iteration for W
#' @param UV_init_sd numeric. SD used to randomly initialize variance of U and V.
#' @param a_u, b_u, a_s, b_s, a_v, b_v, a_r, b_r, a_a, b_a numeric.  Hyper-parameters.
#' @return Object of class 'eClustaLSEC_GEM'.
#' \itemize{
#' \item estimates a list providing the MAP estimates of the model parameters
#' \item algoParms a list providing the optimization parameters
#' \item userInputs a list providing the number of clusters, dimension of the latent space, and the hyperparameters supplied by the user
#' }
#' @import Rcpp
#' @import igraph
#' @import cluster
#' @import rARPACK
eClustaLSEC_GEM = function(A, K=NULL, p=NULL,
                           maxIter=1e3, maxIterVB=100,
                           eps=1e-4, QNSteps=25,
                           CGSteps=25,
                           UV_init_sd=1,
                           a_u = 1,b_u = 1,a_s=1, b_s=1,
                           a_v = a_u,b_v = b_u,a_r=a_s,b_r=b_s,
                           a_a = 1, b_a = 200){

  ### Create objects
  M = ecount(A)
  n = vcount(A)
  EE = cbind(ends(A,1:ecount(A),FALSE))

  cat("\n Indexing edgelist \n")
  temp = indexEdges(EE,n)
  Mi1 = temp$Mi1[,1:max(temp$Mi1Index)]
  Mi2 = temp$Mi2[,1:max(temp$Mi2Index)]
  Mi1Index = temp$Mi1Index
  Mi2Index = temp$Mi2Index
  rm(temp)


  #--- Helper functions:
  {
    .wrappers = list()
    .wrappers$UV = function(x){
      QUV(SS = x[1:n],
          RR = x[n+n*p + 1:n],
          U=matrix(x[n + 1:(n*p)],n,p),
          V=matrix(x[2*n+n*p + 1:(n*p)],n,p),
          W,
          Pmk,Pk,Pmki1,Pmki2,
          tauS, tauR, tauU, tauV,
          EE)
    }
    .wrappers$UVGrad = function(x){
      dQUV(SS = x[1:n],
           RR = x[n+n*p + 1:n],
           U=matrix(x[n + 1:(n*p)],n,p),
           V=matrix(x[2*n+n*p + 1:(n*p)],n,p),
           W,
           Pmk,Pk,Pmki1,Pmki2,
           tauS, tauR, tauU, tauV,
           EE,Mi1Index,Mi1)
    }
    .wrappers$W = function(x){
      QW(S,R,U,V,W=matrix(x,K,p),
         Pmk=Pmk,Pk,
         EE)
    }
    .wrappers$WGrad = function(x){
      c(dQW(S,R,U,V,W=matrix(x,K,p),
            Pmk,Pk,
            EE))
    }
  }
  #--- End helper functions

  ### Initialize
  cat("\n Initialization \n")
  tauV = 1/UV_init_sd^2
  tauU = 1/UV_init_sd^2
  tauR = tauS = 1

  S = rnorm(n); R = rnorm(n)
  S = S - mean(S)
  R = R - mean(R)

  AM <- get.adjacency(A)
  AM_svd <- rARPACK::svds(AM, p)
  U <- AM_svd$u %*% diag(sqrt(AM_svd$d))
  V <- AM_svd$v %*% diag(sqrt(AM_svd$d))

  W = matrix(rnorm(K*p), K, p)

  alpha = a_a/b_a

  oldZ <- rep(1, M)

  cat("\n Beginning VBEM algorithm \n")

  pb = txtProgressBar(0,maxIter-1,style=3)
  for(it in 2:maxIter){

    ### Update Pmk, Pk and alphaTilde (E-step)

    # Initialize VB
    Pmk = computePmk(S,R,U,V,W,rep(1,K),EE)

    # 1st iteration VB
    alphaTld = alpha + colSums(Pmk)
    Elog_pi_k = digamma(alphaTld) - digamma(K*alpha + M)
    Pmk = computePmk(S,R,U,V,W,exp(Elog_pi_k),EE)

    ELBO_Estep = numeric(maxIterVB)
    ELBO_Estep[1] = computeELBO(S,R,U,V,W,Pmk,colSums(Pmk),
                                Elog_pi_k,alphaTld,alpha,EE)

    for(itVB in 2:maxIterVB){
      # Update alphaTilde
      alphaTld_old = alphaTld
      alphaTld = alpha + colSums(Pmk)

      # Update Pmk
      Pmk_old = Pmk
      Elog_pi_k = digamma(alphaTld) - digamma(K*alpha + M)
      Pmk = computePmk(S,R,U,V,W,exp(Elog_pi_k),EE)

      # Check for convergence
      ELBO_Estep[itVB] = computeELBO(S,R,U,V,W,Pmk,colSums(Pmk),
                                     Elog_pi_k,alphaTld,alpha,EE)
      if( ((ELBO_Estep[itVB] - ELBO_Estep[itVB-1])/abs(ELBO_Estep[itVB-1]) < eps ) &
          ( ELBO_Estep[itVB] > ELBO_Estep[itVB-1] ) ) break
    }

    Pk = colSums(Pmk)
    Pmki1 = getPmki(Pmk,Mi1Index,Mi1)
    Pmki2 = getPmki(Pmk,Mi2Index,Mi2)

    newZ <- apply(Pmk, 1, which.max)

    ##Check for convergence
    if(sum(oldZ !=  newZ) < (M/1e3)) break
    oldZ <- newZ

    ### Update U and V
    Opt = optim(par=c(S,c(U),R,c(V)),
                fn = .wrappers$UV,
                gr = .wrappers$UVGrad,
                method="CG",
                control = list(maxit=CGSteps,fnscale=-1))
    S = Opt$par[1:n]
    U = matrix(Opt$par[n + 1:(n*p)],n,p)
    R = Opt$par[n + n*p + 1:n]
    V = matrix(Opt$par[2*n +n*p + 1:(n*p)],n,p)

    S = S - mean(S)
    R = R - mean(R)

    ### Update tauS, tauR, tauU, tauV

    tauS = (a_s + n - 2)/(b_s + sum(S^2))
    tauR = (a_r + n - 2)/(b_r + sum(R^2))
    tauU = (a_u + n*p - 2)/(b_u + sum(U^2))
    tauV = (a_v + n*p - 2)/(b_v + sum(V^2))


    ### Update W
    Opt = optim(par = c(W),
                fn = .wrappers$W,
                gr = .wrappers$WGrad,
                method = "BFGS",
                control = list(maxit=QNSteps,fnscale=-1))
    W = matrix(Opt$par,K,p)

    ### Update alpha
    alpha <- exp(optimize(function(l){
      a <- exp(l)
      obj <- (a_a - 1)*log(a) - b_a*a + lgamma(K*a) - K*lgamma(a) +
        sum((a + Pk - 1)*Elog_pi_k)
      return(obj)
    }, c(-10,10), maximum = T)$maximum)

    setTxtProgressBar(pb,it)
  }


  # Transform the parameter estimates appropriately
  S = S - mean(S)
  R = R - mean(R)
  W_sd = sd(c(W))
  W = W/W_sd

  U = U*W_sd
  V = V*W_sd

  tauU = tauU/W_sd^2
  tauV = tauV/W_sd^2

  ret = list(estimates = list(S=S,
                              R=R,
                              U = U,
                              V = V,
                              W = W,
                              Pmk = Pmk,
                              tauU = tauU,
                              tauV = tauV,
                              tauS = tauS,
                              tauR = tauR,
                              alpha = alpha),
             algoParms = list(eps = eps,
                              numIter = it,
                              maxIter = maxIter),
             userInputs = list(K = K,
                               p = p,
                               a_u = a_u,
                               b_u = b_u,
                               a_v = a_v,
                               b_v = b_v,
                               a_s = a_s,
                               b_s = b_s,
                               a_r = a_r,
                               b_r = b_r,
                               a_a = a_a,
                               b_a = b_a))
  class(ret) = "eClustaLSEC_GEM"
  return(ret)
}
