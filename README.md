## Overview

This package provides functions to implement
a variational Bayes generalized expectation-maximization strategy and 
a Hamiltonian Monte Carlo-within Gibbs algorithm
to estimate the automated latent space 
edge clustering LSEC (aLSEC) model efficiently.

## Installation

``` r
devtools::install_github('hanhtdpham/aLSEC')
```

## Example fitting aLSEC model for UK Faculty network

``` r
# Preamble ----------------------------------------------------------------
library(igraph)
library(scatterplot3d)
library(aLSEC)


# Data --------------------------------------------------------------------
data(UKfaculty,package="igraphdata")
A = UKfaculty
A = remove.edge.attribute(A,"weight")

# Fitting aLSEC model with GEM --------------------------------------------
# Use 15 random starting points
n_starting_points = 15
temp_fits = list()
for(i in 1:n_starting_points){
  set.seed(i)
  temp_fits[[i]] = eClustaLSEC_GEM(A,K=10,p=3)
}

# Choose the best model using ICL
ICLs <- sapply(temp_fits, function(fit) {ICL.eClustaLSEC_GEM(fit, A)})
fit_GEM <- temp_fits[[which.max(ICLs)]]

# Obtain the cluster assignment results
Zhat_GEM <- apply(fit_GEM$estimates$Pmk, 1, which.max)
Zhat_GEM <- as.numeric(factor(Zhat_GEM))


# Visualizing results of GEM ----------------------------------------------

### Code from igraph package to add vertex shape ###
# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)
## end ###

set.seed(1)
lo <- layout_nicely(A)

plot(A,
     layout = lo,
     vertex.shape = c("triangle", "circle", "square", "circle")[V(A)$Group],
     vertex.color=c("white", "white", "white", "black")[V(A)$Group],
     vertex.frame.color = "black",
     vertex.label = NA,
     vertex.size = 3,
     edge.color = RColorBrewer::brewer.pal(4, "Set1")[Zhat_GEM],
     edge.arrow.size = 0)



# Fitting aLSEC model with HMC-within Gibbs -------------------------------

# Draw 30,000 MCMC samples and post process
set.seed(123)
fit_HMC <- eClustaLSEC_HMC(A, GEMObj = NULL,
                           K=10,p=3,
                           stepSize = NULL,
                           numSteps = 25,
                           nSims=3e4,
                           Thin = 5,
                           showPB = TRUE,
                           tuningMethod = "auto",
                           autoTune_n = 250,
                           autoTune_bounds = 0.65 + c(-1,1)*0.1,
                           max_autoTune_count = 100)

fit_HMC_pp <- eClustaLSEC_PostProcess(fit_HMC,
                                      burnin=10000,
                                      thin=1,
                                      summary_function=mean,
                                      A = A,
                                      showPB=TRUE)

# Visualizing MCMC results ------------------------------------------------

# Posterior dist of K (after 10k burn-in but before post-processing)
K_post <- apply(fit_HMC$n_k[(1e4+1):3e4,], 1, function(x) sum(x>0))
barplot(table(K_post), cex.axis=1.5, cex.names=1.5, xlab="K")

# MAP latent feature U and cluster assignment
Uhat  = fit_HMC_pp$HMC_obj$U[which.max(fit_HMC_pp$logPosterior),,]
Zhat  = fit_HMC_pp$HMC_obj$z[which.max(fit_HMC_pp$logPosterior),]
s3d <- scatterplot3d(Uhat,
                     tick.marks=FALSE,
                     xlab='',ylab='',zlab='',
                     pch=c(2:0, 19)[V(A)$Group],
                     cex.symbols=1.5,box=F,
                     col.axis=gray(0.4),
                     mar=rep(0,4),
                     angle=145)
orig  <- s3d$xyz.convert(Uhat)
eList <- as_edgelist(A)
eCols4 <- RColorBrewer::brewer.pal(4, "Set1")[Zhat]
segments(x0 = orig$x[eList[,1]],
         y0 = orig$y[eList[,1]],
         x1 = orig$x[eList[,2]],
         y1 = orig$y[eList[,2]],
         col=adjustcolor(eCols4, alpha.f=0.5))
s3d$points3d(Uhat[V(A)$Group==1,],pch=2,cex=1.5)
s3d$points3d(Uhat[V(A)$Group==2,],pch=1,cex=1.5)
s3d$points3d(Uhat[V(A)$Group==3,],pch=0,cex=1.5)
s3d$points3d(Uhat[V(A)$Group==4,],pch=19,cex=1.5)


# Heat map showing probability of two edges belonging to the same clusters
# Permute edges connecting ppl in the same group to be near each other (block)
# Exclude individuals with unknown schools (actors belonging to Group 4)
eList  <- cbind(ends(A,1:ecount(A),FALSE))
eListG <- eList
for(i in 1:nrow(eList)){
  for(j in 1:2){
    eListG[i,j] <- V(A)$Group[V(A) == eListG[i,j]]
  }
}

Zblock <- c(which(eListG[,1] == 1 & eListG[,1] == eListG[,2]),
            which(eListG[,1] == 2 & eListG[,1] == eListG[,2]),
            which(eListG[,1] == 3 & eListG[,1] == eListG[,2]),
            which(eListG[,1] != eListG[,2] & eListG[,1] != 4 & eListG[,2] != 4))
Z_permuted <- fit_HMC_pp$HMC_obj$z[,Zblock]

nsamp <- dim(fit_HMC_pp$HMC_obj$z)[1]
nedge <- dim(Z_permuted)[2]
P <- diag(1, nedge)
for(i in 1:(nedge - 1)){
  for(j in (i+1):nedge){
    P[i,j] <- sum(Z_permuted[,i] == Z_permuted[,j])/nsamp
    P[j,i] <- P[i,j]
  }
}

gplots::heatmap.2(x=P, Colv = F, Rowv = F,  symm =T, dendrogram = "none", 
                  trace="none", density.info = "none",
                  col = heat.colors(999,rev=T)[-c(1:100)],
                  labRow=F, labCol=F, 
                  margins=c(3,0),
                  keysize=1,
                  key.par=list(mar=c(3.5,0,3,0)),
                  lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, 5), lwid=c(1, 10, 1))                                       
```
