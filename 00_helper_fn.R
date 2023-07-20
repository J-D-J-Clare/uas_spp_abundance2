# === scaling by 2 SD function
scalefn <- function(x) {(x - mean(x))/(2 * sd(x))}

# === notin operator
`%notin%` <- negate(`%in%`)

# === MAE & MAPE functions
mae <- function(x, y) { 
  if(length(x) != length(y)) { stop2("Lengths of vectors differ") }
  return( sum(abs(x - y))/length(x) )
  }
mape <- function(x, y) { 
  if(length(x) != length(y)) { stop2("Lengths of vectors differ") }
  return( sum(abs( (x - y)/y ))/length(x) )
}
smape <- function(x, y) { 
  if(length(x) != length(y)) { stop2("Lengths of vectors differ") }
  return( sum(abs( (x - y)/( (x + y)/2 ) ))/length(x) )
}

# === prepare data inputs for exaxt-sparse-car model in nimble
# Max Joseph post: https://github.com/mbjoseph/CARstan
# brms source code:https://github.com/paul-buerkner/brms/blob/a70f9760b8a5a611dec2ea37cb21227106bb8204/R/data-predictor.R
# This is adapted from brms source code: 
# get_car_terms <- function(W = NULL, gr = NULL, grNV = NULL) {
#   if(length(rownames(W))==0) {
#     stop2("W does not have row names")
#   }
#   # info for car
#   M <- W
#   M <- Matrix::Matrix(M, sparse = TRUE)
#   colnames(M) <- rownames(M) 
#   
#   locations <- levels(factor(gr))
#   Nloc <- length(locations)
#   Jloc <- as.array(match(gr, locations))
#   JlocNV <- as.array(match(grNV, locations))
#   
#   M <- M[locations, locations, drop = FALSE]
#   # plot(rast(as.matrix(M)))
#   
#   edges_rows <- (Matrix::tril(M)@i + 1)
#   edges_cols <- sort(Matrix::triu(M)@i + 1) ## sort to make consistent with rows
#   edges <- cbind("rows" = edges_rows, "cols" = edges_cols)
#   
#   out <- list(
#     Nloc = Nloc, Jloc = Jloc, JlocNV = JlocNV, 
#     Nedges = length(edges_rows),
#     edges1 = as.array(edges_rows),
#     edges2 = as.array(edges_cols)
#   )
#   # additional info for escar
#   Nneigh <- Matrix::colSums(M)
#   
#   inv_sqrt_D <- diag(1 / sqrt(Nneigh))
#   eigenMcar <- t(inv_sqrt_D) %*% M %*% inv_sqrt_D
#   eigenMcar <- eigen(eigenMcar, TRUE, only.values = TRUE)$values
#   
#   out$Nneigh <- Nneigh 
#   out$eigenMcar <- eigenMcar
#   
#   return(out)
# }


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
# desparse_car_lpdf <- nimbleFunction(
#   run = function(x = double(1), car = double(0), 
#                  sdcar = double(0), Nloc = integer(0),
#                  Nedges = integer(0), Nneigh = double(1), 
#                  eigenW = double(1),
#                  edges1 = double(1), edges2 = double(1), log = integer(0, default = 0)) {
#     returnType(double(0))
#     
#     tau <- 1/sdcar^2
#     phit_D <- x * Nneigh
#     phit_W <- nimRep(0, Nloc)
#     ll <- rep(0, Nloc)
#     for (i in 1:Nedges) {
#       phit_W[edges1[i]] <- phit_W[edges1[i]] + x[edges2[i]]
#       phit_W[edges2[i]] <- phit_W[edges2[i]] + x[edges1[i]]
#     }
#     for (i in 1:Nloc) {
#       ll[i] = 0.5 * (log(tau) + log(1 - car * eigenW[i]) -
#                        tau * (phit_D[i] * x[i] - car * (phit_W[i] * x[i]))) 
#     }
# 
#     if(log) return(sum(ll))
#     else return(exp(sum(ll))) 
#   }
# )
# 
# dicar_lpdf <- nimbleFunction(
#   run = function(x = double(1), N = integer(0), 
#                  node1 = double(1), node2 = double(1), log = integer(0, default = 0)) {
#     returnType(double(0))
# 
#     diff <- x[node1] - x[node2]
#     ll <- -0.5 * inprod(diff, diff) + dnorm(sum(x), 0, sd = 0.001 * N, log = TRUE)
# 
#     if(log) return(ll)
#     else return(exp(ll)) 
#   }
# )

# === posterior predictions
predFP <- function(post = NULL, data = NULL, counts = TRUE, val = TRUE) {
  post |> 
    select(contains('beta') & contains('1]'), contains('beta0FP')) -> pars
  if(val == TRUE) { dat <- Data$XVal |> as.matrix() } 
  else { dat <- bind_rows(data.frame(data$XNaive, idx = data$ucellidxNV), 
                     data.frame(data$XVal, idx = data$ucellidxV)) |> 
    arrange(idx) |> as.matrix()}
  
  if(counts == TRUE) {
    preds <- sapply(1:nrow(post), function(i) { 
      rpois(nrow(dat), exp(dat[,1:3] %*% t(pars[i,1:3]) + pars[i, 4] ) ) }) # as.numeric(pars[i, data$siteV + 3])
  } else {
    preds <- sapply(1:nrow(post), function(i) { 
      exp(dat[,1:3] %*% t(pars[i,1:3]) + pars[i, 4] ) }) # as.numeric(pars[i, data$sitefull + 3])
  }
  
  return(preds)
}

predUnk <- function(post = NULL, data = NULL) {
  df <- Data$XNaive
  n <- nrow(df)
  
  # --- TP
  post |> 
    select(contains('beta') & contains('2]'), -contains("beta0FP"), contains("theta"), "sigma") |> as.matrix() -> pars
  post |> 
    select(contains('phi')) |> select(data$ucellidxNV) |> as.matrix() -> re
  post |> 
    select( starts_with("p["), -contains("beta")) |> as.matrix() -> ps
  # --- FP
  post |> 
    select(contains('beta') & contains('1]'), "beta0FP") |> as.matrix() -> parsFP
  # post |> 
  #   select(contains('beta0FP[')) |> as.matrix() -> reFP
  # reFP <- reFP[, data$siteN]

  muFP <- sapply(1:nrow(post), function(i) { exp(df[,1:3] %*% parsFP[i,1:3] + parsFP[i,4]) }) |> t()
  TPhat <- sapply(1:nrow(post), function(i) { 
      mu <- exp( df[,1:3] %*% pars[i,1:3] + pars[i,4]*exp(-df[,4]^2*pars[i,5]) + re[i,] ) 
      ct <- rnbinom(n, prob = mu/(pars[i,6]*mu^2 + mu), size = mu^2/((pars[i,6]*mu^2 + mu) - mu))
    }) |> t()
  
  preds <- sapply(1:n, function(i) {
    rpois(nrow(post), muFP[,i]) + rbinom(nrow(post), TPhat[,i], ps[,data$siteN[i]]) })
  
  return(preds)
}

predTP <- function(post = NULL, data = NULL, counts = TRUE, val = TRUE) {
  post |> 
    select(contains('beta'), contains("theta"), -contains("0FP"), "sigma") |> 
    select(4:9) |> as.matrix() -> pars
  re <- post |> select(contains("phi")) |> as.matrix()
  
  if(val == TRUE) { 
    dat <- data$XVal |> as.matrix() 
    lam <- t(exp(sapply(1:nrow(pars), function(i) { 
      dat[,1:3] %*% pars[i,1:3] + 
        pars[i,4]*exp(-dat[,4]^2*pars[i,5]) + 
        re[i,data$ucellidxV] })) )
    } else { 
      dat <- bind_rows(data.frame(data$XNaive, idx = data$ucellidxNV), 
                          data.frame(data$XVal, idx = data$ucellidxV)) |> 
        arrange(idx) |> as.matrix()
      lam <- t(exp(sapply(1:nrow(pars), function(i) { 
        dat[,1:3] %*% pars[i,1:3] + 
          pars[i,4]*exp(-dat[,4]^2*pars[i,5]) + 
          re[i,] })) )
    }
  preds <- t(sapply(1:nrow(lam), function(i) { 
    s <- pars[i,6]
    return(rnbinom(ncol(lam), prob = lam[i,]/(s* lam[i,]^2 + lam[i,]), size = lam[i,]^2/((lam[i,] + s*lam[i,]^2) - lam[i,]))) }) )
  
  colnames(preds) <- 1:ncol(lam); colnames(lam) <- 1:ncol(lam)
  
  if(counts == TRUE) {
    return(preds)  
  } else {
    return(lam)
  }
}

