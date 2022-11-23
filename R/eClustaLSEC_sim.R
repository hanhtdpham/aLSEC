#' Simulate a network according to the latent space model
#'
#' @param p integer. Number of latent dimensions
#' @param K integer. True number of clusters
#' @param n integer. Number of actors in the network
#' @param M integer. Number of edges in the network
#' @param conc concentration parameter of the mixtures of von Mises-Fisher distributions to generate U, V
#' @param UVShape shape parameter of the Gamma distribution to generate the magnitude of U, V
#' @param UVRate rate parameter of the Gamma distribution to generate the magnitude of U, V
#'
#' @import movMF
#' @import mvtnorm
eClustaLSEC_Sim = function(p = 2,
                           K = 2*p,
                           n = 400,
                           M = round(n * (n-1) * 0.05),
                           conc = c(1.5,2.5,6,25,15,40,100,
                                    1.75,5,10,50,10,15,30)[7*(p-2) + K - 1],
                           UVShape=c(20, 60)[(K>5)+1],
                           UVRate=4){

  library(igraph)

  tauSR = 0.5
  Alpha = rep(1/K,K)

  if(p == 2){
    phi = seq(0,2*pi,l = K + 1)[-(K+1)]
    W = cbind(cos(phi),sin(phi))
    rm(phi)
  }else{
    # Create function to generate points roughly equidistant on sphere
    fibonacci = function(N){
      pts = matrix(0.0, N, 3)
      phi = pi * (3 - sqrt(5))

      for(i in 0:(N-1)){
        y = 1 - (i / (N-1)) * 2
        radius = sqrt(1 - y^2)
        theta = phi * i

        x = cos(theta) * radius
        z = sin(theta) * radius

        pts[i+1,] = c(x,y,z)
      }

      pts
    }

    W = fibonacci(K)
  }

  UDir = W[sample(K,n,TRUE),]
  U = V = matrix(0,n,p)
  UVMagnitude = rgamma(n,shape=UVShape,rate=UVRate)
  for(i in 1:n){
    UVTheta = movMF::rmovMF(2,theta=conc*UDir[i,])
    U[i,] = UVMagnitude[i] * UVTheta[1,]
    V[i,] = UVMagnitude[i] * UVTheta[2,]
  }
  SR = mvtnorm::rmvnorm(n,sigma=tauSR*outer(1:2,1:2,function(x,y)0.75^abs(x-y)))
  eUiWk = exp(SR[,1]%*%matrix(1,1,K) + tcrossprod(U,W))
  eViWk = exp(SR[,2]%*%matrix(1,1,K) + tcrossprod(V,W))

  EL = matrix(NA,10*M,2,
              dimnames = list(NULL,c("sender","receiver")))
  Z = sample(K,size=10*M,replace=T,prob=Alpha)
  nVec = 1:n
  for(m in 1:(10*M)){
    EL[m,1] = sample(n,1,prob=eUiWk[,Z[m]])
    EL[m,2] = sample(nVec[-EL[m,1]],1,prob=eViWk[-EL[m,1],Z[m]])
  }

  edge_distn =
    dplyr::count(as.data.frame(EL),sender,receiver)
  index_to_keep = sample(1:nrow(edge_distn),
                         min(M,nrow(edge_distn)),
                         FALSE,
                         edge_distn$n)

  edge_distn = edge_distn[index_to_keep, ]
  edge_distn$sr = paste(edge_distn$sender, ",", edge_distn$receiver)

  Z = Z[!duplicated(EL)]
  EL = EL[!duplicated(EL),]

  temp = data.frame("sender" = EL[,"sender"], "receiver" = EL[,"receiver"],
                    "Z" = Z)
  temp$sr = paste(temp$sender, ",", temp$receiver)
  temp = temp[temp$sr %in% edge_distn$sr, ]

  EL = as.matrix(temp[,c("sender", "receiver")])
  Z = temp[,"Z", drop=T]
  A = graph_from_edgelist(EL,directed=TRUE)
  E(A)$cluster = Z

  A = delete_vertices(A, which(degree(A) == 0))
  gc()
  return(list(A = A, S = SR[,1], R = SR[,2], U = U, V = V, W = W))
}
