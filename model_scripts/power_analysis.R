pkgs <- c("dplyr", "nimble","rsample")
sapply(pkgs, require, character.only = TRUE)
# source("00_helper_fn.R")
# ---


# This script runs the power analysis.
# The models below are the iterations of the full model, with increasing 
# proportions of validated cells designated as unvalidated. 

# proportions: [0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]

# The following datasets will cary the indices:
# 0. generate the split/subset of thedata where val == 1
# 1. nTP, FN, FP, Xval, siteV, ucellidxV [1:729]
# 2. n, XNaive, siteN, ucellidxNV [1:5876]

set.seed(12345)


for(i in 10:10){
  datin <- sort(rep(c("", "025"), 10))
  idx <- rep(1:10, 2)  

  props <- c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
  # 2. load site level data
  datsite <- read.csv("detection/uas_spp_abundance2/data/celldat_site.csv")
  
  # load abundance
  dat <- read.csv(paste0("detection/uas_spp_abundance2/data/celldat_input", datin[i], ".csv")) 
  # ================================================================
  
  split <- lapply(props, function(x) { 
    training(initial_split( filter(dat, val == 1), prop = x, strata = site)) |> pull(ucid) })
  
  if(i %in% c(1, 10)) { dat <- dat}
	else { dat <- dat |>  mutate(val = ifelse(ucid %in% split[[idx[i]]], 1, 0)) }
  
  # ================================================================
  
  # === icar input
  carinfo <- readRDS("detection/uas_spp_abundance2/data/icar_data_8.rds") 
  
  # adj <- unlist(carinfo$nblist)
  adj <- sapply(list.files("detection/uas_spp_abundance2/data/icar_data_site/", pattern = "_8_", full.names = TRUE), function(x){ unlist(readRDS(x)$nblist) }) |> unlist()
  weights <- adj/adj
  num <- sapply(list.files("detection/uas_spp_abundance2/data/icar_data_site/", pattern = "_8_", full.names = TRUE), function(x){ lengths(readRDS(x)$nblist) }) |> unlist()
  L <- length(adj)
  # --- add site index for phi[]
  siteadj <- sapply(list.files("detection/uas_spp_abundance2/data/icar_data_site/", pattern = "_8_", full.names = TRUE), function(x){ length( unlist(readRDS(x)$nblist) ) })
  sitephi <- table(dat$site)
  
  idxsiteadj <- matrix(rep(0, 20), nr = 10)
  idxsiteadj[,1] <- cumsum(siteadj) - siteadj + 1
  idxsiteadj[,2] <- cumsum(siteadj)
  
  idxsitephi <- matrix(rep(0, 20), nr = 10)
  idxsitephi[,1] <- cumsum(sitephi) - sitephi + 1
  idxsitephi[,2] <- cumsum(sitephi)
  
  # ===
  
  ucellidxV <- dat |>  
    mutate(uidx = 1:n()) |> filter(val == 1) |> pull(uidx) # , site == !!site
  ucellidxNV <- dat |> 
    mutate(uidx = 1:n()) |> filter(val == 0) |> pull(uidx) # , site == !!site
  # ======
  
  # ====== 
  set.seed(12354)
  Data <- list(
    nTP = dat |> mutate(nTP = TP + FN) |> filter(val == 1) |> pull(nTP),
    FN = dat |> filter(val == 1) |> pull(FN),
    FP = dat |> filter(val == 1) |> pull(FP),
    n = dat |> filter(val == 0) |> pull(Unk),
    
    XNaive = dat |> filter(val == 0) |> 
      mutate(Int = 1) |>
      select(tri_sc, tpi_sc, aspect_sc, dist) |> as.matrix(),
    XVal = dat |> filter(val == 1) |> 
      mutate(Int = 1) |>
      select(tri_sc, tpi_sc, aspect_sc, dist) |> as.matrix(), 
    Xsite = datsite |>
      mutate(Int = 1) |>
      select(Int, wind_max_sc, avght_sc, seg_sc) |> as.matrix(),
    bunrt = dat |> pull(burnt),
    
    siteV = dat |> filter(val == 1) |> pull(site) |> as.factor() |> as.numeric(),
    siteN = dat |> filter(val == 0) |> pull(site) |> as.factor() |> as.numeric(),
    sitefull = dat |> pull(site) |> as.factor() |> as.numeric(),
    ucellidxV = ucellidxV,
    ucellidxNV = ucellidxNV,
    # --- 
    
    Nloc = nrow(dat), 
    Nedges = carinfo$N_edges,
    node1 = carinfo$node1, 
    node2 = carinfo$node2,
    adj = adj,
    num = num,
    weights = weights, 
    L = L,
    idxsiteadj = idxsiteadj,
    idxsitephi = idxsitephi
  )
  # ============================================================================
run_MCMC_allcode <- function(seed, Data = NULL) {
  library(nimble)
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
    })
  assign('dCompoundNMix', dCompoundNMix, envir = .GlobalEnv)
  assign('rCompoundNMix', rCompoundNMix, envir = .GlobalEnv)

  # ---
  Mod <- nimbleCode(
    {
      for (b in 1:2) {
        for (j in 1:nxf) {
          beta[j, b] ~ dnorm(0, sd = .5)
        }
      }
      
      for(b in 1:nsf){
        psi[b] ~ dnorm(0, sd = .5)
      }
      
      beta0FP ~ dnorm(0, sd = 1)
      theta1 ~ T(dnorm(0, sd = 1), 0, )
      theta2 ~ dgamma(3, 0.58)
      sigma ~ T(dnorm(0, 5), 0, )
      
      
      for(s in 1:nsite){
        logit(p[s]) <- inprod(psi[1:nsf], Xsite[s, 1:nsf])
        carsd[s] ~ T(dnorm(0, sd = 1), 0, )
      }
      
      for (j in 1:nNaive) {
        log(mu[j, 1]) <- inprod(beta[1:nxf, 1], XNaive[j, 1:nxf]) + 
          beta0FP
        log(mu[j, 2]) <- inprod(beta[1:nxf, 2], XNaive[j, 1:nxf]) + 
          theta1 * exp(-XNaive[j, 4] * theta2) + phi[ucellidxNV[j]]
        
        NTP[j] ~ dnegbin(prob = mu[j, 2]/(mu[j, 2] + sigma*mu[j, 2]^2), size = mu[j, 2]^2/((mu[j, 2] + sigma*mu[j, 2]^2) - mu[j, 2])) 
        n[j] ~ dCompoundNMix(N = NTP[j], p = p[sidx2[j]], lamFP = mu[j, 1])
      }
      
      for (j in 1:nVal) {
        log(mu2[j, 1]) <- inprod(beta[1:nxf, 1], XVal[j, 1:nxf]) + 
          beta0FP
        log(mu2[j, 2]) <- inprod(beta[1:nxf, 2], XVal[j, 1:nxf]) + 
          theta1 * exp(-XVal[j, 4] * theta2) + phi[ucellidxV[j]]
        
        FP[j] ~ dpois(mu2[j, 1])
        nTP[j] ~ dnegbin(prob = mu2[j, 2]/(mu2[j, 2] + sigma*mu2[j, 2]^2), size = mu2[j, 2]^2/((mu2[j, 2] + sigma*mu2[j, 2]^2) - mu2[j, 2])) 
        FN[j] ~ dbin(1 - p[sidx[j]], nTP[j])
      }
      
      for(i in 1:nsite){
        phi[idxsitephi[i,1]:idxsitephi[i,2]] ~ dcar_normal(adj = adj[idxsiteadj[i,1]:idxsiteadj[i,2]], weights = weights[idxsiteadj[i,1]:idxsiteadj[i,2]], tau = carsd[i], num = num[idxsitephi[i,1]:idxsitephi[i,2]], zero_mean = 0)
      }
    }
  )
  
  Constants <- list(  nxf  = 4 - 1, 
                      nNaive = nrow(Data$XNaive), 
                      nVal = nrow(Data$XVal), 
                      nsite = nrow(Data$Xsite),
                      sidx = Data$siteV,
                      sidx2 = Data$siteN,
                      nsf = 4,
                      # ---
                      Nloc = Data$Nloc,
                      ucellidxV = Data$ucellidxV, 
                      ucellidxNV = Data$ucellidxNV,
                      L = Data$L,
                      # --- 
                      idxsiteadj = Data$idxsiteadj,
                      idxsitephi = Data$idxsitephi
  )
  
  Inits <- list(
    beta0FP = 0,
    beta = matrix(0, nr = Constants$nxf, nc = 2),
    NTP = round(Data$n/2), 
    phi = rep(0, Constants$Nloc),
    theta1 = 1,
    theta2 = 1,
    sigma = 1,
    carsd = rep(.1, Constants$nsite),
    psi = rep(0, Constants$nsf)
  )
  
  
  Shrub<- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                      data=Data, inits = Inits)
  
  
  ##monitors -- what things to trace
  ShrubConf <- configureMCMC(Shrub, monitors = c("beta0FP", "beta", "sigma", "p", "psi", "phi", "carsd", "theta1", "theta2", "NTP", "mu", "mu2"), 
                            onlySlice = FALSE)
  
  
  Rmcmc <- buildMCMC(ShrubConf)
  compMCMC <- compileNimble(Rmcmc, Shrub)
  
  samp <- runMCMC(mcmc = compMCMC$Rmcmc,
                  niter = 220000, nburnin = 20000, thin = 200, 
                  nchains = 1, samplesAsCodaMCMC = TRUE)
  
  return(samp)
  }

  library(parallel)
  cl <- makeCluster(8)
  lsamp <- parLapply(cl = cl, X = 1:8, fun = run_MCMC_allcode, Data = Data)
  stopCluster(cl)

  if(i <= 10) {j = 1} else {j = 2}
  saveRDS(list(posterior = lsamp), 
           paste0("detection/models/power_analysis/nf", j,".2.1nb_power_", rep(props, 2)[i], ".rds") )
  # saveRDS(list(data = Data, df = dat), 
  #	paste0("detection/models/power_analysis/data", j, ".2.1nb_power_data_", rep(props, 2)[i], ".rds"))
  
}

#stop cluster


