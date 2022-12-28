#' Automated edge clustering via Hamiltonian Monte Carlo
#'
#' 'eClustaLSEC_HMC' implements Hamiltonian Monte Carlo (HMC) for clustering edges of a network using a latent space approach
#'
#' @param A igraph object
#' @param GEMObj Object of class 'eClustaLSEC_GEM'.  If not supplied, the GEM algorithm
#' will run to initialize the HMC.
#' @param K integer. Upper truncation on the number of clusters.  Will be extracted from GEMObj if supplied.
#' @param p integer. Dimension of latent space.  Will be extracted from GEMObj if supplied.
#' @param stepSize numeric. Step size for HMC.
#' @param numSteps integer. Number of leapfrog steps for HMC.
#' @param nSims integer. Number of posterior draws to be obtained through HMC
#' @param Thin integer. Only each \code{Thin}-th sample will be saved.
#' @param showPB logical. Whether or not to show the progress bar.
#' @param tuningMethod One of 'auto' or 'none'. If 'auto', pilot runs will be used to
#'  obtain an acceptance rate within \code{autoTune_bounds}.
#' @param autoTune_n integer. The number of posterior samples used in each pilot run to
#' estimate the acceptance rate if auto-tuning.
#' @param autoTune_bounds numeric vector. The lower and upper bounds of desired acceptance
#' rates if auto-tuning.  Note that an acceptance rate during a pilot run does not guarantee
#' an overall acceptance rate.  This is NOT an adaptive MCMC algorithm.
#' @param max_autoTune_count integer.  The maximum number of pilot runs allowed if auto-tuning.
#' @param a_u, b_u, a_s, b_s, a_v, b_v, a_r, b_r, a_a, b_a numeric.  Hyperparameters.
#' @return Object of class 'eClustaLSEC_HMC'.
#'
eClustaLSEC_HMC = function(A,GEMObj= NULL,
                           K=NULL,p=NULL,
                           stepSize = NULL,
                           numSteps = 25,
                           nSims=1e3,
                           Thin = 5,
                           showPB = TRUE,
                           tuningMethod = c("auto","none"),
                           autoTune_n = 250,
                           autoTune_bounds = 0.65 + c(-1,1)*0.1,
                           max_autoTune_count = 100,
                           a_u = 0.1,b_u = 0.1,a_v=a_u, b_v=b_u,
                           a_s = 0.1, b_s = 0.1, a_r = a_s, b_r = b_s,
                           a_a = 1, b_a = 200){

  ## Start of HMC function ----------------------------------------------
  tuningMethod = tuningMethod[1]
  if(is.null(GEMObj)){
    GEMObj =
      eClustaLSEC_GEM(A,K=K,p=p,maxIter = 100, maxIterVB=100,
                      eps=1e-4, QNSteps=25, CGSteps = 25,
                      UV_init_sd=1,
                      # alpha_init=NULL,
                      a_u = a_u,b_u = b_u, a_v = a_v, b_v=b_v,
                      a_s = a_s, b_s = b_s, a_r=a_r, b_r=b_r,
                      a_a = a_a, b_a = b_a)
  }

  ### Create objects
  M = ecount(A)
  n = vcount(A)
  p = ncol(GEMObj$estimates$U)
  K = nrow(GEMObj$estimates$W)
  EE = cbind(ends(A,1:ecount(A),FALSE))

  SS = RR = matrix(0.0,nSims,n)
  UU = VV = array(0.0,c(nSims,n,p))
  WW = array(0.0,c(nSims,K,p))
  tauS = tauR = tauU = tauV = lambda = numeric(nSims)
  n_k = matrix(0.0,nSims,K)
  pi = matrix(0.0,nSims,K)
  z = matrix(0L,nSims,M)

  # cat("\n Indexing edgelist \n")
  temp = indexEdges(EE,n)
  Mi1 = temp$Mi1[,1:max(temp$Mi1Index)]
  Mi2 = temp$Mi2[,1:max(temp$Mi2Index)]
  Mi1Index = temp$Mi1Index
  Mi2Index = temp$Mi2Index
  rm(temp)

  ### Initialize
  z[1,] = apply(GEMObj$estimates$Pmk,1,which.max)
  n_k[1, ] = sapply(1:K,function(k) sum(z[1,] == k))
  SS[1,] = GEMObj$estimates$S
  RR[1,] = GEMObj$estimates$R
  UU[1,,] = GEMObj$estimates$U
  VV[1,,] = GEMObj$estimates$V
  WW[1,,] = GEMObj$estimates$W
  pi[1,] = GEMObj$estimates$alphaTld/sum(GEMObj$estimates$alphaTld)
  lambda[1] = log(GEMObj$estimates$alpha)
  tauS[1] = (GEMObj$userInputs$a_s + n)/(GEMObj$userInputs$b_s + sum(SS[1,]^2))
  tauR[1] = (GEMObj$userInputs$a_r + n)/(GEMObj$userInputs$b_r + sum(RR[1,]^2))
  tauU[1] = GEMObj$estimates$tauU
  tauV[1] = GEMObj$estimates$tauV

  if(is.null(stepSize)){
    stepSize = 0.1/(n*(2*p + 2) + K*p)^(1/3) #Roberts and Rosenthal (1998)
  }


  drawUVWl = function(iter = 1,ind = 1){
    # See, e.g., Ch5 in MCMC Handbook, Fig 5.2

    # Helper function
    vec_grad = function(obj){
      c(obj$UV[,1],obj$UV[,p+2],c(obj$UV[,1 + 1:p]),c(obj$UV[,p + 2 + 1:p]),c(obj$W), obj$lambda)
    }


    S_new = SS[iter,]
    R_new = RR[iter,]
    U_new = UU[iter,,]
    V_new = VV[iter,,]
    W_new = WW[iter,,]
    lambda_new = lambda[iter]

    momentum_1 <- momentum_0 <- rnorm(n*(2*p + 2) + K*p + 1)

    # Make a half step for momentum at the beginning
    grad_0 =
      dPostUVW(S_new,R_new,U_new,V_new,W_new,z[iter,],
               tauS[iter],tauR[iter],tauU[iter],tauV[iter],EE,n_k[iter,],
               Mi1Index,Mi1,Mi2Index,Mi2)
    grad_0[["lambda"]] = exp(lambda_new)*sum(digamma(n_k[iter,] + exp(lambda_new))) +
      K*exp(lambda_new)*(digamma(K*exp(lambda_new)) - digamma(exp(lambda_new)) -
                           digamma(K*exp(lambda_new) + M)) +
      a_a - b_a*exp(lambda_new)

    momentum_1 = momentum_1 +
      0.5*stepSize[ind]*vec_grad(grad_0)

    # Alternate full steps for position and momentum
    for(ell in 1:numSteps){
      # Make a full step for the position
      S_new = S_new +
        stepSize[ind]*momentum_1[1:n]
      R_new = R_new +
        stepSize[ind]*momentum_1[n + 1:n]
      U_new = U_new +
        stepSize[ind]*matrix(momentum_1[2*n + 1:(n*p)],n,p)
      V_new = V_new +
        stepSize[ind]*matrix(momentum_1[2*n + n*p + 1:(n*p)],n,p)
      W_new = W_new +
        stepSize[ind]*matrix(momentum_1[2*n + 2*n*p + 1:(K*p)],K,p)
      lambda_new = lambda_new +
        stepSize[ind]*momentum_1[2*n + 2*n*p + K*p + 1]

      grad_0 =
        dPostUVW(S_new,R_new,U_new,V_new,W_new,z[iter,],
                 tauS[iter],tauR[iter],tauU[iter],tauV[iter],EE,n_k[iter,] ,
                 Mi1Index,Mi1,Mi2Index,Mi2)
      grad_0[["lambda"]] =exp(lambda_new)*sum(digamma(n_k[iter,] + exp(lambda_new))) +
        K*exp(lambda_new)*(digamma(K*exp(lambda_new)) - digamma(exp(lambda_new)) -
                             digamma(K*exp(lambda_new) + M)) +
        a_a - b_a*exp(lambda_new)

      # Make a full step for the momentum, except at end of trajectory
      if(ell != numSteps){
        momentum_1 = momentum_1 + stepSize[ind]*vec_grad(grad_0)
      }

    }

    # Make a half step for momentum at the end.
    momentum_1 = momentum_1 +
      0.5*stepSize[ind]*vec_grad(grad_0)

    # Negate momentum at end of trajectory to make the proposal symmetric
    momentum_1 = -momentum_1

    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    logAccProb =
      evalConditionalPost(z[iter,],S_new,R_new,U_new,V_new,W_new,
                          tauS[iter],tauR[iter],tauU[iter],tauV[iter],EE,
                          n_k[iter,] , lambda_new, a_a, b_a) -
      evalConditionalPost(z[iter,],SS[iter,],RR[iter,],UU[iter,,],VV[iter,,],WW[iter,,],
                          tauS[iter],tauR[iter],tauU[iter],tauV[iter],EE,
                          n_k[iter,] , lambda[iter], a_a, b_a) -
      0.5*(sum(momentum_1^2) - sum(momentum_0^2))

    if(runif(1) < exp(logAccProb)){
      return(list(S=S_new, R=R_new, U=U_new, V=V_new, W=W_new, lambda=lambda_new, accept=1))
    }else{
      return(list(S=SS[iter,], R=RR[iter,], U=UU[iter,,], V=VV[iter,,], W=WW[iter,,], lambda=lambda[iter],accept=0))
    }
  }

  #################
  ### Auto-Tune ###
  #################
  if(tolower(tuningMethod) == "auto"){
    cat("-------------------------\n")
    cat("--- Begin auto-tuning ---\n")
    cat("-------------------------\n")
    accRate = 0
    autoTune_count = 0
    stepSize = c(stepSize,NA,NA)
    turnNumber = 1
    gr = (1+sqrt(5))*0.5
    while( ((accRate < autoTune_bounds[1]) | (accRate > autoTune_bounds[2])) & (autoTune_count < max_autoTune_count)){

      ### First turn of auto-tuning
      if(turnNumber == 1){
        accRate = 0
        # if(showPB) pb = txtProgressBar(0,autoTune_n,style=3)
        for(it in 1:autoTune_n){
          new_draws = drawUVWl(iter=1,ind=1)
          accRate = accRate + new_draws$accept/autoTune_n
          if(new_draws$accept > 0 ){
            SS[1,] = new_draws$S
            RR[1,] = new_draws$R
            UU[1,,] = new_draws$U
            VV[1,,] = new_draws$V
            WW[1,,] = new_draws$W
            lambda[1] = new_draws$lambda
          }
          # if(showPB) setTxtProgressBar(pb,it)
        }
        cat(paste0("--- Acceptance Rate = ",round(accRate,4)," ---\n"))
        autoTune_count = 1
      }

      ### Second turn of auto-tuning
      if(turnNumber == 2){
        stepSize[2] = stepSize[1]
        if(accRate < autoTune_bounds[1]){
          stepSize[2] = stepSize[2]*0.5
          while( (accRate < autoTune_bounds[1]) & (autoTune_count < max_autoTune_count)){
            accRate = 0
            # if(showPB) pb = txtProgressBar(0,autoTune_n,style=3)
            for(it in 1:autoTune_n){
              new_draws = drawUVWl(iter=1,ind=2)
              accRate = accRate + new_draws$accept/autoTune_n
              if(new_draws$accept > 0 ){
                SS[1,] = new_draws$S
                RR[1,] = new_draws$R
                UU[1,,] = new_draws$U
                VV[1,,] = new_draws$V
                WW[1,,] = new_draws$W
                lambda[1] = new_draws$lambda
              }
              # if(showPB) setTxtProgressBar(pb,it)
            }
            cat(paste0("--- Acceptance Rate = ",round(accRate,4)," ---\n"))
            autoTune_count = autoTune_count + 1
            if(accRate < autoTune_bounds[1]) stepSize[2] = stepSize[2]*0.5
          }
        }else{
          stepSize[2] = stepSize[2]*2
          while( (accRate > autoTune_bounds[2]) & (autoTune_count < max_autoTune_count)){
            accRate = 0
            # if(showPB) pb = txtProgressBar(0,autoTune_n,style=3)
            for(it in 1:autoTune_n){
              new_draws = drawUVWl(iter=1,ind=2)
              accRate = accRate + new_draws$accept/autoTune_n
              if(new_draws$accept > 0 ){
                SS[1,] = new_draws$S
                RR[1,] = new_draws$R
                UU[1,,] = new_draws$U
                VV[1,,] = new_draws$V
                WW[1,,] = new_draws$W
                lambda[1] = new_draws$lambda
              }
              # if(showPB) setTxtProgressBar(pb,it)
            }
            cat(paste0("--- Acceptance Rate = ",round(accRate,4)," ---\n"))
            autoTune_count = autoTune_count + 1
            if(accRate > autoTune_bounds[2]) stepSize[2] = stepSize[2]*2
          }
        }
      }

      ### Past second turn of auto-tuning
      if(turnNumber > 2){
        stepSize[1:2] = sort(stepSize[1:2])
        stepSize[3] = (stepSize[1] + gr*stepSize[2])/(1 + gr) #Golden ratio search

        accRate = 0
        # if(showPB) pb = txtProgressBar(0,autoTune_n,style=3)
        for(it in 1:autoTune_n){
          new_draws = drawUVWl(iter=1,ind=3)
          accRate = accRate + new_draws$accept/autoTune_n
          if(new_draws$accept > 0 ){
            SS[1,] = new_draws$S
            RR[1,] = new_draws$R
            UU[1,,] = new_draws$U
            VV[1,,] = new_draws$V
            WW[1,,] = new_draws$W
            lambda[1] = new_draws$lambda
          }
          # if(showPB) setTxtProgressBar(pb,it)
        }
        cat(paste0("--- Acceptance Rate = ",round(accRate,4)," ---\n"))
        autoTune_count = autoTune_count + 1
        if(accRate > autoTune_bounds[2]){
          stepSize[1] = stepSize[3]
        }else{
          stepSize[2] = stepSize[3]
        }
      }

      turnNumber = turnNumber + 1
      if((accRate > autoTune_bounds[1]) & (accRate < autoTune_bounds[2])){
        autoTune_count = Inf
        if(turnNumber == 2) stepSize = stepSize[1]
        if(turnNumber == 3) stepSize = stepSize[2]
        if(turnNumber > 3) stepSize = stepSize[3]
      }
    }#End: while
  }


  ###############
  ### Run HMC ###
  ###############
  accRate = 0
  if(showPB) pb = txtProgressBar(0,nSims*Thin,style=3)
  for(it in 2:nSims){
    z[it,] = z[it-1,]
    SS[it, ] = SS[it-1,]
    RR[it,] = RR[it-1,]
    UU[it,,] = UU[it-1,,]
    VV[it,,] = VV[it-1,,]
    WW[it,,] = WW[it-1,,]
    pi[it,]  = pi[it-1,]
    n_k[it,] = n_k[it-1,]
    lambda[it] = lambda[it-1]
    tauS[it] = tauS[it-1]
    tauR[it] = tauR[it-1]
    tauU[it] = tauU[it-1]
    tauV[it] = tauV[it-1]

    for(th in 1:Thin){

      { #Draw U,V, and W, and lambda
        new_draws = drawUVWl(iter=it,ind=1)
        accRate = accRate + new_draws$accept/Thin/nSims
        if(new_draws$accept > 0 ){
          SS[it,] = new_draws$S
          RR[it,] = new_draws$R
          UU[it,,] = new_draws$U
          VV[it,,] = new_draws$V
          WW[it,,] = new_draws$W
          lambda[it] = new_draws$lambda
        }
      }
      { #Draw z and pi's
        pi[it,] = rdirichlet(1, exp(lambda[it]) + n_k[it, ])
        z[it,] = drawZ(SS[it,],RR[it,],UU[it,,],VV[it,,],WW[it,,],pi[it,],EE)
        n_k[it, ] = sapply(1:K,function(k) sum(z[it,] == k))
      }

      { #Draw the tau's
        tauS[it] = rgamma(1,
                          shape = 0.5*(GEMObj$userInputs$a_s + n),
                          rate = 0.5*(GEMObj$userInputs$b_s + sum(SS[it,]^2)))
        tauR[it] = rgamma(1,
                          shape = 0.5*(GEMObj$userInputs$a_r + n),
                          rate = 0.5*(GEMObj$userInputs$b_r + sum(RR[it,]^2)))
        tauU[it] = rgamma(1,
                          shape = 0.5*(GEMObj$userInputs$a_u + n*p),
                          rate = 0.5*(GEMObj$userInputs$b_u + sum(UU[it,,]^2)))
        tauV[it] = rgamma(1,
                          shape = 0.5*(GEMObj$userInputs$a_v + n*p),
                          rate = 0.5*(GEMObj$userInputs$b_v + sum(VV[it,,]^2)))
      }

      if(showPB) setTxtProgressBar(pb,(it-1)*Thin + th)
    }
  }

  ret = list(z=z,S=SS,R=RR,U=UU,V=VV,W=WW,n_k=n_k,pi=pi,alpha=exp(lambda),
             tauS=tauS,tauR=tauR,tauU=tauU,tauV=tauV,
             K=K, nSims = nSims, p=p, thinning=Thin,
             acceptanceRate = accRate, stepSize=stepSize,
             GEMObj = GEMObj)

  class(ret) = "eClustaLSEC_HMC"
  return(ret)
}
