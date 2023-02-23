pkgs <- c("tidyverse", "sf", "rstan", "terra", 'nimble')
sapply(pkgs, require, character.only = TRUE)
scalefn <- function(x) {(x - mean(x))/(2 * sd(x))}
# ---

# === data is currently set up to attempt capture-recapture
# n_obs: number of drone detected plants 
# ngt25: number of field mapped plants taller than 25 cm
# n: number of field mapped plants
# nTP: number of field mapped + detected with UAS
# nTPM: number of field mapped + detected + not NTP b/c 1 segment = 2 field points
# nFP: number of field mapped + detected + classified as ARTR but not in fact ARTR
# can be assumed to be absent, since we have only 17 plants across all sites
# nM: number of field plants not detected with UAS
# nU: number of plants detected w/ UAS, but unknown validation status


# load cell-level data
dat0 <- read.csv("data/celldat_full.csv")  |>
  mutate(dist_sq = dist^2) |>
  filter(!is.na(tpi) | !is.na(slope))  # |> filter(val == 1) |> select(-val)

# ---
site <- unique(dat0$site)[1]
# ---
# === covariates
# het3: structural heterogeneity at coarse scale
# het6: structural heterogeneity at medium scale
# het8: structural heterogeneity at fine scale
# tpi: topographic position index
# slope: slope 
# elev: cell elevation
# max.ht: maximum canopy height in a cell
# avg.ht: average canopy height
# dist: distance to fire line (negative = inside unburnt) 
# dist_sq: squared distance - likely makes sense

# --- add standardized predictors
dat0 |> 
  filter(site == !!site) |>
  mutate(across(c(het3, het6, het8, tpi, slope, elev, max.ht, avg.ht, dist, dist_sq), 
                scalefn, .names = "{.col}_sc")) -> dat

# --- corplot for nTP: True Positives counts
cr <- dat |> filter(val == 1) |> 
  select(where(is.numeric), -contains("_sc"), -val) |> cor()
corrplot::corrplot(cr)

# --- 
# --- mixture likelihood function
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
# ---

Mod<-nimbleCode({
    for (b in 1:2){ 
      for(j in 1:nxf){ # # of covariates
        beta[j, b] ~ dnorm(0, sd = 2) ###describe intensity of point processes
      }
    }
  p ~ dbeta(2, 2) ###pr real shrub detected from air
  
  for (j in 1:length(nNaive)) {
    log(mu[j, 1]) <- inprod(beta[1:nxf, 1], XNaive[j, 1:nxf] )  # FP
    log(mu[j, 2]) <- inprod(beta[1:nxf, 2], XNaive[j, 1:nxf] )  # TP
    NTP[j]~dpois(mu[j, 2])
    n[j]~dCompoundNMix(N=NTP[j], p=p, lamFP=mu[j, 1])
  }
  
  for (j in 1:length(nVal)) {
    mu2[j, 1] <- exp(inprod(beta[1:nxf,1], XVal[j,1:nxf]))
    nTP[j]~dpois(mu2[j, 2])
    mu2[j, 2] <- exp(inprod(beta[1:nxf,2], XVal[j,1:nxf]))
    nFP[j]~dpois(mu2[j, 1])
  }
  
  ###validated plants...
  ###guess at the number of unvalided shrubs that exist
  EN<-sum(mu[1:nNaive, 2])+sum(mu2[1:nVal, 2])
  N<-sum(NTP[1:nNaive])+sum(nTP[1:nVal])
})

Data <- list(
  n = dat |> filter(val == 0) |> pull(n_obs),
  nTP = dat |> filter(val == 1) |> pull(ngt25),
  nFP = dat |> filter(val == 1) |> 
    mutate(ns = n_obs - ngt25, 
           ns = ifelse(ns < 0, 0, ns)) |> pull(ns),
  
  XNaive= dat |> filter(val == 0) |> 
    mutate(Int = 1) |>
    select(Int, burn, het3_sc, het8_sc, max.ht_sc, dist_sc) |> as.matrix(),
  XVal= dat |> filter(val == 1) |> 
    mutate(Int = 1) |>
    select(Int, burn, het3_sc, het8_sc, max.ht_sc, dist_sc) |> as.matrix()
  )

Constants <- list(  nxf  = 6, 
                    nNaive = nrow(dat) - sum(dat$val), 
                    nVal = sum(dat$val) )

Inits <- list(p = .5, 
              beta = matrix(0, nr = Constants$nxf, nc = 2),
              NTP = round(Data$n*.5), 
              mu = matrix(1, nr = Constants$nNaive, nc = 2) ,
              mu2 = matrix(1, nr = Constants$nVal, nc = 2) 
              )


Shrub<- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                    data=Data, inits = Inits)

##monitors--what things to trace
ShrubConf<- configureMCMC(Shrub, monitors = c("beta", "N", "EN"))


Rmcmc <- buildMCMC(ShrubConf)
compMCMC <- compileNimble(Rmcmc, Shrub)

samps <- runMCMC(mcmc = compMCMC$Rmcmc,
               niter = 50000, nburnin = 40000, thin = 10, 
               nchains = 1, samplesAsCodaMCMC = TRUE)

pred <- rpois(Constants$nVal,
              exp( Data$XVal %*% apply(samps[,1:Constants$nxf], 2, mean) ) + 
                exp(Data$XVal %*% apply(samps[,(Constants$nxf+1):(2*ncol(Data$XVal))], 2, mean) ) ) 
pred <- exp(Data$XVal %*% apply(samps[,1:Constants$nxf], 2, mean)) + 
  exp(Data$XVal %*% apply(samps[,(Constants$nxf+1):(2*ncol(Data$XVal))], 2, mean) )
plot(log(pred) ~ log(Data$n), pch = 19, cex = .5)
abline(0, 1)
cor(pred, Data$nTP)^2

dat |> 
  filter(val == 0) |>
  mutate(pred = #rpois(Constants$nNaive, 
                      exp(Data$XNaive[,1:2] %*% apply(samps[,3:4], 2, mean)) ) |>
  bind_rows(dat |> 
              filter(val == 1) |>
              mutate(pred = #rpois(Constants$nVal, 
                                  exp(Data$XVal[,1:2] %*% apply(samps[,1:2], 2, mean)) )) |>
  st_as_sf(coords = c("east", "north"), crs = 32611) |>
  # mutate(nTP = log(nTP)) |>
  filter(site == !!site) -> out

plot(out[,'pred'], pch = 19)
plot(out[,'nTP'], pch = 19)

