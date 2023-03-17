pkgs <- c("tidyverse", "sf", "rstan", "terra", 'nimble')
sapply(pkgs, require, character.only = TRUE)
scalefn <- function(x) {(x - mean(x))/(2 * sd(x))}
`%notin%` <- negate(`%in%`)
# ---

# === data is currently set up to attempt capture-recapture
# note: only > 25 cm plants considered & classified as ARTR from random-forest
# see power point for 'rules' of counts
# TP: count of field mapped + matched with UAS (verified cells only)
# FP: count of UAS polygons not matched with field pts (verified cells only)
# FN: count of field mapped plant + not detected with UAS (verified cells only)
# Unk: count of plants detected w/ UAS (unvalidated cells only)

# load covariates
datcov <- read.csv("data/celldat_covar.csv") |> 
  mutate(dist_sq = dist^2) |>
  filter(!is.na(tpi) | !is.na(slope)) 

# load counts
datct <- read.csv("data/celldat_match_counts.csv") |> 
  filter(ucid %in% datcov$ucid)

dat0 <- datcov |> 
  bind_cols(datct |> select(TP, FN, FP, Unk))

# ---
site <- unique(dat0$site)[3]
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
  filter(site == !!site, burnt == 0) |> #, 
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
    
    # n[j]~dpois(p*mu[j, 2] + mu[j, 1])
    n[j]~dCompoundNMix(N=NTP[j], p=p, lamFP=mu[j, 1])
  }
  
  for (j in 1:length(nVal)) {
    log(mu2[j, 1]) <- inprod(beta[1:nxf,1], XVal[j,1:nxf]) # FP
    nFP[j]~dpois(mu2[j, 1])
    log(mu2[j, 2]) <- inprod(beta[1:nxf,2], XVal[j,1:nxf]) # TP
    nTP[j]~dpois(mu2[j, 2])
    FN[j]~binom(1-p, nTP[j])
  }
  
  ###validated plants...
  ###guess at the number of unvalided shrubs that exist
  EN<-sum(mu[1:nNaive, 2])+sum(mu2[1:nVal, 2])
  N<-sum(NTP[1:nNaive])+sum(nTP[1:nVal])
  NVal <- sum(mu2[1:nVal, 2])
})

Data <- list(
  nTP = dat |> mutate(nTP = TP + FN) |> filter(val == 1) |> pull(nTP),
  nFP = dat |> filter(val == 1) |> pull(FP),
  n = dat |> filter(val == 0) |> pull(Unk),
  
  XNaive= dat |> filter(val == 0) |> 
    mutate(Int = 1) |>
    select(Int, max.ht_sc) |> as.matrix(),
  XVal= dat |> filter(val == 1) |> 
    mutate(Int = 1) |>
    select(Int, max.ht_sc) |> as.matrix()
  )

Constants <- list(  nxf  = 2, 
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
ShrubConf<- configureMCMC(Shrub, monitors = c("beta", "N", "EN", "NVal", "p"), 
                          onlySlice = FALSE)


Rmcmc <- buildMCMC(ShrubConf)
compMCMC <- compileNimble(Rmcmc, Shrub)

samp <- runMCMC(mcmc = compMCMC$Rmcmc,
               niter = 100000, nburnin = 50000, thin = 100, 
               nchains = 4, samplesAsCodaMCMC = TRUE)
# ---
library(MCMCvis)
MCMCtrace(samp, pdf = FALSE, params = "beta")
MCMCtrace(samp, pdf = FALSE, params = "p")
MCMCtrace(samp, pdf = FALSE, params = c("EN", "N","NVal") )
cat("nTP: ", sum(Data$nTP))
cat("nUnk: ", sum(Data$n))
# ---

samps <- map(samp, as.data.frame) |> bind_rows()
pred <- exp(Data$XVal %*% apply(samps[,6:7], 2, mean))
plot(pred ~ (Data$nTP), pch = 19, cex = .5)
abline(0, 1)
cor(pred, Data$nTP)^2

dat |>
  filter(val == 0) |>
  mutate(pred = exp(Data$XNaive %*% apply(samps[,6:7], 2, mean)) ) |>
  bind_rows(dat |> 
              filter(val == 1) |>
              mutate(pred = exp(Data$XVal %*% apply(samps[,6:7], 2, mean)) )) |>
  st_as_sf(coords = c("east", "north"), crs = 32611) |>
  filter(site == !!site) -> out

plot(stars::st_as_stars(out[,'pred']), pch = 19)
plot(out[,'TP'], pch = 19)

# === 
pred <- exp(apply(samps[,6:7], 1, function(x) { Data$XVal %*% x } ) )
pred |> 
  as.data.frame() |> 
  mutate(nTP = Data$nTP, burnt = Data$XVal[,2], id = 1:nrow(Data$XVal) ) |>
  pivot_longer(cols = contains("V")) |> select(-name) |>
  group_by(id) |> summarize_all(median) |> ungroup() |>
  ggplot(aes(x = nTP, y = value, colour = burnt)) + geom_point(alpha = .5)  
  # geom_abline(intercept = 0, slope = 1)


