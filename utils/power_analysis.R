pkgs <- c("tidyverse", "sf", "rsample")
sapply(pkgs, require, character.only = TRUE)
source("00_helper_fn.R")
# ---


# This script runs the power analysis.
# The models below are the iterations of the full model, with increasing 
# proportions of validated cells designated as unvalidated. 

# proportions: [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]

# The following datasets will cary the indices:
# 0. generate the split/subset of thedata where val == 1
# 1. nTP, FN, FP, Xval, siteV, ucellidxV [1:729]
# 2. n, XNaive, siteN, ucellidxNV [1:5876]

set.seed(12345)

# 1. 
props <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)


for(i in 1:length(props)){
  # 2. load site level data
  datsite <- read.csv("data/site_level_data.csv") |> 
    mutate(seg = pareaseg/denseg) |>
    select(site, elev, slope, jdate, wind_max, wind_avg, gps_err, gcp_err,
           pareasegartr, densegartr, pareaseg, denseg, seg, avght) |> 
    mutate(across(c(elev, slope, jdate, wind_max, wind_avg, gps_err, gcp_err,
                    pareasegartr, densegartr, pareaseg, denseg, seg, avght), scalefn, 
                  .names = "{.col}_sc"))
  
  # load covariates
  datcov <- read.csv("data/celldat_covar.csv") |> 
    filter(!is.na(tpi) | !is.na(slope)) |> 
    mutate(dist = ifelse(dist < 0, 0, dist/100 ) ) |>
    select(-c(het3, het6, het8, max.ht, avg.ht))  |>
    mutate(across(c(aspect, aspect.sm), 
                  function(x){(180 - (x - 225))/360}))
  
  # load counts
  datct <- read.csv("data/celldat_match_counts.csv") |> 
    filter(ucid %in% datcov$ucid)
  
  dat0 <- datcov |> 
    bind_cols(datct |> select(TP, FN, FP, Unk)) |> 
    mutate(nTP = TP + FN, dist = dist)
  
  # ================================================================
  
  split <- lapply(props, function(x) { 
    training(initial_split( filter(dat0, val == 1), prop = x, strata = site)) |> pull(ucid) })
  
  dat0 |> 
    mutate(val = ifelse(ucid %in% split[[i]], 1, 0)) -> dat0
  
  # ================================================================
  
  # --- add standardized predictors
  dat0 |> 
    as.data.frame() |>
    mutate(across(c(tri, tri.sm, tpi, tpi.sm, slope, slope.sm, aspect, aspect.sm, dist), 
                  scalefn, .names = "{.col}_sc")) -> dat
  
  # === icar input
  carinfo <- readRDS("data/icar_data_8.rds") 
  
  adj <- unlist(carinfo$nblist)
  weights <- adj/adj
  num <- lengths(carinfo$nblist)
  L <- length(adj)
  
  # --- add site index for phi[]
  siteadj <- sapply(list.files("data/icar_data_site/", pattern = "_8_", full.names = TRUE), function(x){ length( unlist(readRDS(x)$nblist) ) })
  sitephi <- table(dat$site)
  
  idxsiteadj <- matrix(rep(0, 20), nr = 10)
  idxsiteadj[,1] <- cumsum(siteadj) - siteadj + 1
  idxsiteadj[,2] <- cumsum(siteadj)
  
  idxsitephi <- matrix(rep(0, 20), nr = 10)
  idxsitephi[,1] <- cumsum(sitephi) - sitephi + 1
  idxsitephi[,2] <- cumsum(sitephi)
  
  # ===
  ucellidxV <- dat0 |> filter(site %in% unique(dat$site)) |> 
    mutate(uidx = 1:n()) |> filter(val == 1) |> pull(uidx) 
  ucellidxNV <- dat0 |> filter(site %in% unique(dat$site)) |> 
    mutate(uidx = 1:n()) |> filter(val == 0) |> pull(uidx) 
  # ===================================================================
  
  
  # ====== 
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
      beta0 ~ dnorm(0, sd = 1)
      theta1 ~ dexp(2.5)
      theta2 ~ dgamma(3, 0.58)
      sigma_inv ~ dnorm(0, 5)
      sigma <- 1/sigma_inv
      
      
      for(s in 1:nsite){
        logit(p[s]) <- inprod(psi[1:nsf], Xsite[s, 1:nsf])
        carsd[s] ~ T(dnorm(0, sd = 1), 0, )
      }
      
      for (j in 1:nNaive) {
        log(mu[j, 1]) <- inprod(beta[1:nxf, 1], XNaive[j, 1:nxf]) + 
          beta0FP
        log(mu[j, 2]) <- inprod(beta[1:nxf, 2], XNaive[j, 1:nxf]) + 
          theta1 * exp(-XNaive[j, 4]^2 * theta2) + phi[ucellidxNV[j]]
        
        NTP[j] ~ dnegbin(prob = sigma/(mu[j, 2] + sigma), size = sigma)
        n[j] ~ dCompoundNMix(N = NTP[j], p = p[sidx2[j]], lamFP = mu[j, 1])
      }
      
      for (j in 1:nVal) {
        log(mu2[j, 1]) <- inprod(beta[1:nxf, 1], XVal[j, 1:nxf]) + 
          beta0FP
        FP[j] ~ dpois(mu2[j, 1])
        log(mu2[j, 2]) <- inprod(beta[1:nxf, 2], XVal[j, 1:nxf]) + 
          theta1 * exp(-XVal[j, 4]^2 * theta2) + phi[ucellidxV[j]]
        
        nTP[j] ~ dnegbin(prob = sigma/(mu[j, 2] + sigma), size = sigma)
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
    beta0 = 0,
    beta0FP = 0,
    beta = matrix(0, nr = Constants$nxf, nc = 2),
    NTP = round(Data$n/2), 
    phi = rep(0, Constants$Nloc),
    theta1 = 1,
    theta2 = 1,
    sigma_inv = 1,
    carsd = rep(.1, Constants$nsite),
    psi = rep(0, Constants$nsf)
  )
  
  
  Shrub<- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                      data=Data, inits = Inits)
  
  
  ##monitors--what things to trace
  ShrubConf<- configureMCMC(Shrub, monitors = c("beta0FP", "beta", "sigma", "sigma_inv", "p", "psi", "phi", "carsd", "theta1", "theta2", "NTP", "mu", "mu2"), 
                            onlySlice = FALSE)
  
  
  Rmcmc <- buildMCMC(ShrubConf)
  compMCMC <- compileNimble(Rmcmc, Shrub)
  
  samp <- runMCMC(mcmc = compMCMC$Rmcmc,
                  niter = 70000, nburnin = 10000, thin = 30, 
                  nchains = 4, samplesAsCodaMCMC = TRUE)
  
  saveRDS(list(Mod = Mod, Data = Data, Constants = Constants,
               Inits = Inits, nimbleMod = Shrub, ModConf = ShrubConf,
               mcmc = Rmcmc, compiledModel = compMCMC, posterior = samp), 
          paste0("../../models/nimblefit1.2dd_p", props[i], ".rds") )
}



