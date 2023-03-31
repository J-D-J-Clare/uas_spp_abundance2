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
  filter(!is.na(tpi) | !is.na(slope)) |> 
  select(-c(het3, het6, het8, max.ht, avg.ht))

datcov <- st_read("data/celldat_covar.geojson") |> 
  filter(!is.na(tpi) | !is.na(slope)) |> 
  # mutate(dist = ifelse(dist < 0, 0, log(abs(dist+1)) ) ) |>
  mutate(dist = ifelse(dist <= 0, -log(abs(dist-1)), log(abs(dist+1)) ) ) |>
  select(-c(het3, het6, het8, max.ht, avg.ht))  

plot(datcov[datcov$site==site,'burnt'])

# load counts
datct <- read.csv("data/celldat_match_counts.csv") |> 
  filter(ucid %in% datcov$ucid)

dat0 <- datcov |> 
  bind_cols(datct |> select(TP, FN, FP, Unk)) |> 
  mutate(nTP = TP + FN)

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

# ---
site <- unique(dat0$site)[8]
# --- add standardized predictors
dat0 |> 
  as.data.frame() |>
  filter(site == !!site) |> 
  mutate(across(c(tri, tri.sm, tpi, tpi.sm, slope, slope.sm, flowdir, flowdir.sm, dist), 
                scalefn, .names = "{.col}_sc")) -> dat
plot(datcov[datcov$site==site,'dist'])

# --- corplot for nTP: True Positives counts
cr <- dat |> filter(val == 1) |> 
  as.data.frame() |>
  select(contains("_sc"), -val, nTP) |> cor()
corrplot::corrplot(cr)

# ======
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
        beta[j, b] ~ dnorm(0, sd = .5) ###describe intensity of point processes
      }
    }
  p ~ dbeta(2, 2) ###pr real shrub detected from air
  
  for (j in 1:nNaive) {
    log(mu[j, 1]) <- inprod(beta[1:nxf, 1], XNaive[j, 1:nxf] )  # FP
    log(mu[j, 2]) <- inprod(beta[1:nxf, 2], XNaive[j, 1:nxf] )  # TP
    NTP[j]~dpois(mu[j, 2])

    n[j]~dCompoundNMix(N=NTP[j], p=p, lamFP=mu[j, 1])
  }
  
  for (j in 1:nVal) {
    log(mu2[j, 1]) <- inprod(beta[1:nxf,1], XVal[j,1:nxf]) # FP
    FP[j]~dpois(mu2[j, 1])
    log(mu2[j, 2]) <- inprod(beta[1:nxf,2], XVal[j,1:nxf])  # TP
    nTP[j] ~ dpois(mu2[j, 2])
    FN[j]~dbin(1-p, nTP[j])
  }
})

Data <- list(
  nTP = dat |> mutate(nTP = TP + FN) |> filter(val == 1) |> pull(nTP),
  FN = dat |> filter(val == 1) |> pull(FN),
  FP = dat |> filter(val == 1) |> pull(FP),
  n = dat |> filter(val == 0) |> pull(Unk),
  
  XNaive= dat |> filter(val == 0) |> 
    mutate(Int = 1) |>
    select(Int, dist_sc, tpi_sc, tpi.sm_sc) |> as.matrix(),
  XVal= dat |> filter(val == 1) |> 
    mutate(Int = 1) |>
    select(Int, dist_sc, tpi_sc, tpi.sm_sc) |> as.matrix()
  )

Constants <- list(  nxf  = 4, 
                    nNaive = nrow(dat) - sum(dat$val), 
                    nVal = sum(dat$val) )

Inits <- list(p = .5, 
              beta = matrix(0, nr = Constants$nxf, nc = 2),
              NTP = rep(0, Constants$nNaive)
              )


Shrub<- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                    data=Data, inits = Inits)

##monitors--what things to trace
ShrubConf<- configureMCMC(Shrub, monitors = c("beta", "p"), 
                          onlySlice = FALSE)


Rmcmc <- buildMCMC(ShrubConf)
compMCMC <- compileNimble(Rmcmc, Shrub)

samp <- runMCMC(mcmc = compMCMC$Rmcmc,
               niter = 10000, nburnin = 5000, thin = 10, 
               nchains = 4, samplesAsCodaMCMC = TRUE)
# ---
library(MCMCvis)
MCMCtrace(samp, pdf = FALSE, params = "beta")
MCMCtrace(samp, pdf = FALSE, params = "p")

cat("nTP: ", sum(Data$nTP))
cat("nUnk: ", sum(Data$n))
# ---

samps <- map(samp, as.data.frame) |> bind_rows()
pred <- apply(exp(apply(samps[,5:8], 1, function(x) { 
  Data$XVal %*% x })), 1, mean)
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


