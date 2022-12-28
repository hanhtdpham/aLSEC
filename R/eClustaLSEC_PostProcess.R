#' Post-process Hamiltonian Monte Carlo
#'
#' 'eClustaLSEC_PostProcess' post processes the HMC output for the aLSEC model to
#' make summaries of the posterior draws possible.
#'
#' @param HMC_obj Object of class 'eClustaLSEC_HMC'.
#' @param burnin integer. Number of leading posterior draws to remove.
#' @param thin integer. Only each \code{Thin}-th sample will be processed.
#' @param summary_function function, e.g. mean or median. Function used to summarize
#' the posterior samples.
#' @param A igraph.
#' @param showPB logical. Whether or not to show the progress bar.
#' @return Named List.
#' \itemize{
#' \item HMC_obj object of class 'eClustaLSEC_HMC'
#' \item keep_seq integer vector giving the index of the original HMC draws kept in HMC_obj
#' \item HMC_summary object of class 'eClustaLSEC_HMC_summary'
#' \item logPosterior vector giving the log posterior values at each HMC draw.
#' }
eClustaLSEC_PostProcess = function(HMC_obj,
                                   burnin=1,
                                   thin=1,
                                   summary_function=mean,
                                   A = NULL,
                                   showPB=TRUE){
  nSims = nrow(HMC_obj$z)
  n = dim(HMC_obj$U)[2]
  iter_seq = seq(burnin+1,nSims,by=thin)

  cat("--- Remove burin and thin ---\n")
  for(var_name in c("tauS","tauR","tauU","tauV","alpha")){
    HMC_obj[[var_name]] = HMC_obj[[var_name]][iter_seq]
  }
  for(var_name in c("z","S","R","n_k","pi")){
    HMC_obj[[var_name]] = HMC_obj[[var_name]][iter_seq,]
  }
  for(var_name in c("U","V","W")){
    HMC_obj[[var_name]] = HMC_obj[[var_name]][iter_seq,,]
  }

  HMC_obj$burnin = burnin*HMC_obj$thinning
  HMC_obj$thinning = HMC_obj$thinning*thin
  nSims = length(iter_seq)

  # Center S and R
  cat("--- Center R and S ---\n")

  S_new = HMC_obj$S
  R_new = HMC_obj$R
  for(it in 1:nSims){
    S_new[it,] = S_new[it,] - mean(S_new[it,])
    R_new[it,] = R_new[it,] - mean(R_new[it,])
  }

  # Rescale W and subsequently U and V (and recenter), and finally tauU and tauV
  U_new = HMC_obj$U
  V_new = HMC_obj$V
  W_new = HMC_obj$W
  tauU_new = HMC_obj$tauU
  tauV_new = HMC_obj$tauV
  cat("\n")
  cat("--- Rescale W and subsequently U and V (and recenter) and their tau's ---\n")
  if(showPB) pb = txtProgressBar(0,nSims,style=3)
  for(it in 1:nSims){
    W_sd = sd(c(HMC_obj$W[it,,]))
    W_new[it,,] = W_new[it,,]/W_sd

    U_new[it,,] = U_new[it,,]*W_sd
    V_new[it,,] = V_new[it,,]*W_sd

    tauU_new[it] = tauU_new[it]/W_sd^2
    tauV_new[it] = tauV_new[it]/W_sd^2

    if(showPB) setTxtProgressBar(pb,it)
  }
  W_sd = sd(c(HMC_obj$GEMObj$estimates$W))

  # Remove draws where number of cluster is not posterior mode of K
  cat("\n")
  cat("--- Keep draws with number of clusters == posterior mode of K ---\n")
  K_0 = apply(HMC_obj$n_k, 1, function(x) sum(x>0))
  HMC_obj$K = as.numeric(tail(names(sort(table(K_0))),1))

  S_new = S_new[K_0 == HMC_obj$K,]
  R_new = R_new[K_0 == HMC_obj$K,]
  U_new = U_new[K_0 == HMC_obj$K,,]
  V_new = V_new[K_0 == HMC_obj$K,,]
  W_new = W_new[K_0 == HMC_obj$K,,]
  tauU_new = tauU_new[K_0 == HMC_obj$K]
  tauV_new = tauV_new[K_0 == HMC_obj$K]

  for(var_name in c("tauS","tauR")){
    HMC_obj[[var_name]] = HMC_obj[[var_name]][K_0 == HMC_obj$K]
  }
  for(var_name in c("z","n_k","pi")){
    HMC_obj[[var_name]] = HMC_obj[[var_name]][K_0 == HMC_obj$K,]
  }
  nSims = nrow(S_new)


  # Find posterior mode for reference (for both label switching and rotation of U,V,W)
  # Note: integrate out pi from posterior because log pi_k can be undefined when pi_k==0
  cat(" Finding marginal posterior mode for reference\n")
  logPosterior = numeric(nSims)
  EE = cbind(ends(A,1:ecount(A),FALSE))
  if(showPB) pb = txtProgressBar(0,nSims,style=3)
  for(it in 1:nSims){
    logPosterior[it] =
      evalMargPost(
        HMC_obj$z[it,],
        S_new[it,],R_new[it,],
        U_new[it,,],V_new[it,,],W_new[it,,],
        HMC_obj$tauS[it],HMC_obj$tauR[it],
        tauU_new[it],tauV_new[it],
        HMC_obj$n_k[it,], log(HMC_obj$alpha[it]),
        EE,
        HMC_obj$GEMObj$userInputs$a_s,HMC_obj$GEMObj$userInputs$b_s,
        HMC_obj$GEMObj$userInputs$a_r,HMC_obj$GEMObj$userInputs$b_r,
        HMC_obj$GEMObj$userInputs$a_u,HMC_obj$GEMObj$userInputs$a_v,
        HMC_obj$GEMObj$userInputs$b_u,HMC_obj$GEMObj$userInputs$b_v,
        HMC_obj$GEMObj$userInputs$a_a,HMC_obj$GEMObj$userInputs$b_a)

    if(showPB) setTxtProgressBar(pb,it)
  }



  # Remove empty components in each draw
  cat("\n")
  cat("--- Remove parameters associated with empty components ---\n")
  n_k_new <- matrix(0, nSims, HMC_obj$K)
  pi_new  <- matrix(0, nSims, HMC_obj$K)
  HMC_obj$W <- W_new
  W_new   <- array(0.0, c(nSims, HMC_obj$K, dim(W_new)[3]))
  z_new   <- HMC_obj$z
  if(showPB) pb = txtProgressBar(0,nSims,style=3)
  for(it in 1:nSims){
    non.empty    <- which(HMC_obj$n_k[it,] > 0)
    z_new[it,]   <- as.numeric(factor(HMC_obj$z[it,]))
    n_k_new[it,] <- HMC_obj$n_k[it,][non.empty]
    pi_new[it,]  <- HMC_obj$pi[it,][non.empty]
    W_new[it,,]  <- HMC_obj$W[it,non.empty,]
    if(showPB) setTxtProgressBar(pb,it)
  }

  # Relabel z and then W
  cat("\n")
  cat("--- Run Equivalence Classes Representatives (ECR) ---\n")
  zHat = z_new[which.max(logPosterior),]
  ecr_results = ecr(zpivot = zHat, z = z_new, K = HMC_obj$K)
  if(showPB) pb = txtProgressBar(0,nSims,style=3)
  for(it in 1:nSims){
    z_new[it,]  = match(z_new[it,],ecr_results$permutations[it,])
    pi_new[it,] = pi_new[it, ecr_results$permutations[it,]]
    n_k_new[it,] = n_k_new[it, ecr_results$permutations[it,]]
    if(showPB) setTxtProgressBar(pb,it)
  }
  W_new = permute.mcmc(W_new, ecr_results$permutations)$output



  # Rotate rbind(U,V) (and then W)
  cat("\n")
  cat("--- Rotate U and V and subsequently W ---\n")

  UV_reference = rbind(U_new[which.max(logPosterior),,],
                       V_new[which.max(logPosterior),,])

  if(showPB) pb = txtProgressBar(0,nSims,style=3)
  cat("\n")
  cat(" Perform rotations\n")
  for(it in 1:nSims){
    procr = procrustes(X=rbind(U_new[it,,],V_new[it,,]),
                       Xstar=UV_reference,
                       translation=FALSE,dilation=FALSE)
    U_new[it,,] = procr$X.new[1:n,]
    V_new[it,,] = procr$X.new[n + 1:n,]
    W_new[it,,] = W_new[it,,]%*%procr$R
    if(showPB) setTxtProgressBar(pb,it)
  }

  # Recenter U and V
  for(it in 1:nSims){
    U_new[it,,] = scale(U_new[it,,],scale=F)
    V_new[it,,] = scale(V_new[it,,],scale=F)
  }

  # Output and obtain summary
  HMC_obj$z = z_new
  HMC_obj$S = S_new
  HMC_obj$R = R_new
  HMC_obj$U = U_new
  HMC_obj$V = V_new
  HMC_obj$W = W_new
  HMC_obj$tauU = tauU_new
  HMC_obj$tauV = tauV_new
  HMC_obj$n_k = n_k_new
  HMC_obj$pi  = pi_new

  HMC_summary = list()
  HMC_summary$z = matrix(0,ncol(HMC_obj$z),HMC_obj$K)
  for(i in 1:ncol(HMC_obj$z)){
    for(k in 1:HMC_obj$K){
      HMC_summary$z[i,k] = mean(HMC_obj$z[,i] == k)
    }
  }
  HMC_summary$S = apply(HMC_obj$S[1:nSims,],2,summary_function)
  HMC_summary$R = apply(HMC_obj$R[1:nSims,],2,summary_function)
  HMC_summary$U = apply(HMC_obj$U[1:nSims,,],2:3,summary_function)
  HMC_summary$V = apply(HMC_obj$V[1:nSims,,],2:3,summary_function)
  HMC_summary$W = apply(HMC_obj$W[1:nSims,,],2:3,summary_function)
  HMC_summary$n_k = apply(HMC_obj$n_k[1:nSims,],2,summary_function)
  HMC_summary$tauS = summary_function(HMC_obj$tauS[1:nSims])
  HMC_summary$tauR = summary_function(HMC_obj$tauR[1:nSims])
  HMC_summary$tauU = summary_function(HMC_obj$tauU[1:nSims])
  HMC_summary$tauV = summary_function(HMC_obj$tauV[1:nSims])

  class(HMC_summary) = "eClustaLSEC_HMC_summary"

  return(list(HMC_obj = HMC_obj,
              keep_seq = iter_seq,
              HMC_summary = HMC_summary,
              logPosterior = logPosterior))
}
