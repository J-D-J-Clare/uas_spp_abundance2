pkgs <- c("dplyr", 'nimble')
sapply(pkgs, require, character.only = TRUE)
# source("00_helper_fn.R")
# ---
path <- "detection/uas_spp_abundance2/"

# === 
# note: only > or < 25 cm plants (and segments) considered & classified as ARTR from random-forest
# see power point for 'rules' of counts
# TP: count of field mapped + matched with UAS (verified cells only)
# FP: count of UAS polygons not matched with field pts (verified cells only)
# FN: count of field mapped plant & not detected with UAS (verified cells only)
# Unk: count of plants detected w/ UAS (unvalidated cells only)
# ---


# load site level data
datsite <- read.csv(paste0(path, "/data/celldat_site.csv"))

# load abundance data
dat <- read.csv(paste0(path, "/data/celldat_input025.csv"))

# === icar input
carinfo <- readRDS(paste0(path, "/data/icar_data_8.rds") ) # joint
# adj <- unlist(carinfo$nblist)

# icar input by-site
adj <- sapply(list.files(paste0(path, "/data/icar_data_site/"), pattern = "_8_", full.names = TRUE), function(x){ unlist(readRDS(x)$nblist) }) |> unlist()
weights <- adj/adj
num <- sapply(list.files(paste0(path, "/data/icar_data_site/"), pattern = "_8_", full.names = TRUE), function(x){ lengths(readRDS(x)$nblist) }) |> unlist()
L <- length(adj)
# --- add site index for phi[]
siteadj <- sapply(list.files(paste0(path, "/data/icar_data_site/"), pattern = "_8_", full.names = TRUE), function(x){ length( unlist(readRDS(x)$nblist) ) })
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
  burnt = dat |> pull(burnt),
  
  siteV = dat |> filter(val == 1) |> pull(site) |> as.factor() |> as.numeric(),
  siteN = dat |> filter(val == 0) |> pull(site) |> as.factor() |> as.numeric(),
  sitefull = dat |> pull(site) |> as.factor() |> as.numeric(),
  ucellidxV = ucellidxV,
  ucellidxNV = ucellidxNV,
  # --- 
  
  Nloc = nrow(dat), 
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
  rCompoundNMix <- nimbleFunction(
    #In principle, this should generate a random value
    #sensu rnorm(), but an appropriate rng was not needed for the application
    #(rather, the HPC needed a place holding function to operate). If you need to
    #generate random numbers (to impute data, for example), this should be redefined.
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

# ============================================================================
library(parallel)
cl <- makeCluster(8)
lsamp <- parLapply(cl = cl, X = 1:8, 
                   fun = run_MCMC_allcode, 
                   Data = Data)

# It's good practice to close the cluster when you're done with it.
stopCluster(cl)
# ============================================================================

saveRDS(lsamp, "detection/models/nf2.2.1nb_posterior.rds")

# ---
# library(MCMCvis)
# MCMCtrace(lsamp, pdf = FALSE, params = "beta")
# MCMCtrace(lsamp, pdf = FALSE, params = "theta2")
# MCMCsummary(lsamp, params = c("carsd", "sigma", "sigma_inv", "p", "psi", "beta", "theta1", "theta2"))
# 
# phisumm <- MCMCsummary(lsamp, params = "phi")
# hist(phisumm$Rhat)
# ntpsumm <- MCMCsummary(lsamp, params = "NTP")
# hist(ntpsumm$Rhat)


out <- lapply(lsamp, function(x) { return(x |> as.data.frame(x) |> 
                                            select(-contains("mu"), -contains("NTP")) ) }) |> bind_rows()
saveRDS(out, "detection/models/nf2.2.1nb_pars.rds")

# ========= 
saveRDS(lapply(lsamp, function(x) { return(x |> as.data.frame(x) |> 
                                             select(contains("mu")) ) }) |> bind_rows(), 
        "detection/models/nf2.2.1nb_mu.rds")

saveRDS(lapply(lsamp, function(x) { return(x |> as.data.frame(x) |> 
                                             select(contains("NTP")) ) }) |> bind_rows(), 
        "detection/models/nf2.2.1nb_NTP.rds")
