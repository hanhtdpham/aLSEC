## Overview

This package implments a variational Bayes generalized
expectation-maximization strategy to estimate the automated latent space
edge clustering LSEC (aLSEC) model efficiently.

## Installation

``` r
devtools::install_github('hanhtdpham/aLSEC')
```

## Example fitting aLSEC model for UK Faculty network

``` r
# Preamble ----------------------------------------------------------------
library(igraph)
library(MetBrewer)
library(aLSEC)


# Data --------------------------------------------------------------------
data(UKfaculty,package="igraphdata")
A = UKfaculty
A = remove.edge.attribute(A,"weight")

# Fitting aLSEC model -----------------------------------------------------
# Use 15 random starting points
n_starting_points = 15
temp_fits = list()
for(i in 1:n_starting_points){
  set.seed(i)
  temp_fits[[i]] = eClustaLSEC_GEM(A,K=10,p=3)
}

# Choose the best model using ICL
ICLs <- sapply(temp_fits, function(fit) {ICL.eClustaLSEC_GEM(fit, A)})
fit_uk <- temp_fits[[which.max(ICLs)]]

# Obtain the cluster assignment results
Zhat_uk <- apply(fit_uk$estimates$Pmk, 1, which.max)
Zhat_uk <- as.numeric(factor(Zhat_uk))



# Visualizing results  ----------------------------------------------------

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
     vertex.shape = c("circle", "square", "triangle", "circle")[V(A)$Group],
     vertex.color=c("black", "black", "black", "white")[V(A)$Group],
     vertex.frame.color = "black",
     vertex.label = NA,
     vertex.size = 2,
     edge.color = met.brewer("Austria", 10)[Zhat_uk],
     edge.arrow.size = 0)
```
