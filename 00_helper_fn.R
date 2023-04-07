# === scaling function
scalefn <- function(x) {(x - mean(x))/(2 * sd(x))}


# === notin operator
`%notin%` <- negate(`%in%`)


# === prepare data inputs for exaxt-sparse-car model in nimble
# Max Joseph post: https://github.com/mbjoseph/CARstan
# brms source code:https://github.com/paul-buerkner/brms/blob/a70f9760b8a5a611dec2ea37cb21227106bb8204/R/data-predictor.R
# This is adapted from brms source code: 
get_car_terms <- function(W = NULL, gr = NULL) {
  if(length(rownames(W))==0) {
    stop2("W does not have row names")
  }
  # info for car
  M <- W
  M <- Matrix::Matrix(M, sparse = TRUE)
  colnames(M) <- rownames(M) 

  locations <- gr #levels(factor(gr))
  Nloc <- length(locations)
  Jloc <- as.array(match(gr, locations))

  M <- M[locations, locations, drop = FALSE]
  plot(rast(as.matrix(M)))

  edges_rows <- (Matrix::tril(M)@i + 1)
  edges_cols <- sort(Matrix::triu(M)@i + 1) ## sort to make consistent with rows
  edges <- cbind("rows" = edges_rows, "cols" = edges_cols)
  
  out <- list(
    Nloc = Nloc, Jloc = Jloc, Nedges = length(edges_rows),
    edges1 = as.array(edges_rows),
    edges2 = as.array(edges_cols)
  )
  # additional info for escar
  Nneigh <- Matrix::colSums(M)

  inv_sqrt_D <- diag(1 / sqrt(Nneigh))
  eigenMcar <- t(inv_sqrt_D) %*% M %*% inv_sqrt_D
  eigenMcar <- eigen(eigenMcar, TRUE, only.values = TRUE)$values
  
  out$Nneigh <- Nneigh 
  out$eigenMcar <- eigenMcar

  return(out)
}


# === nimble functionality
library(nimble)
# --- compound mixture
dCompoundNMix<-nimbleFunction(
  run = function(x = integer(0), N = integer(0), p=double(0), lamFP=double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    ll1<-nimRep(0, x+1)
    ll2<-nimRep(0, x+1)
    for (i in 1:(x+1)) {
      ll1[i]<-dbinom(i-1, N, p)
      ll2[i]<-dpois(x-i+1, lamFP)
    }
    logProb<-log(sum(ll1*ll2))
    if(log) return(logProb)
    else return(exp(logProb)) 
  })

rCompoundNMix <- nimbleFunction( # this is a placeholder for parallel implementation
  run = function(n = integer(0), N = integer(0), p=double(0), 
                 lamFP=double(0) ) {
    returnType(double(0))
    return(1)
  }
)

# --- escar following Max Josephs specification: https://github.com/mbjoseph/CARstan
# /* Return the log probability of a proper conditional autoregressive (CAR)
# * prior with a sparse representation for the adjacency matrix
# * Full credit to Max Joseph (https://github.com/mbjoseph/CARstan)
# * Args:
# *   phi: Vector containing the CAR parameters for each location
# *   car: Dependence (usually spatial) parameter for the CAR prior
# *   sdcar: Standard deviation parameter for the CAR prior
# *   Nloc: Number of locations
# *   Nedges: Number of edges (adjacency pairs)
# *   Nneigh: Number of neighbors for each location
# *   eigenW: Eigenvalues of D^(-1/2) * W * D^(-1/2)
# *   edges1, edges2: Sparse representation of adjacency matrix
# * Details:
#   *   D = Diag(Nneigh)
# * Returns:
#   *   Log probability density of CAR prior up to additive constant
# */
desparse_car_lpdf <- nimbleFunction(
  run = function(x = double(0), phi = double(0), car = double(0), 
                 sdcar = double(0), Nloc = integer(0),
                 Nedges = integer(0), Nneigh = integer(0), 
                 eigenW = double(0),
                 edges1 = integer(0), edges2 = integer(0)) {
    returnType(double(0))
    
    tau <- 1/sdcar^2
    phit_D <- t(phi * Nneigh)
    phit_W <- t(rep(0, Nloc))
    for (i in 1:Nedges) {
      phit_W[edges1[i]] <- phit_W[edges1[i]] + phi[edges2[i]];
      phit_W[edges2[i]] <- phit_W[edges2[i]] + phi[edges1[i]];
    }
    for (i in 1:Nloc) {
      ldet[i] = log(1 - car * eigenW[i]);
    }
    return(0.5 * (Nloc * log(tau) + sum(ldet) -
                    tau * (phit_D %*% phi - car * (phit_W %*% phi))) )
  }
) 
