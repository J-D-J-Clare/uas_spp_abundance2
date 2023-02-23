###load data
CellDat<-read.csv("/Users/johnclare/Documents/Drone/celldat.csv")
SpDat<-read.csv("/Users/johnclare/Documents/Drone/spdat_joint0.csv")
SpDat$z[SpDat$z==101]<-NA


###a couple of things to consider related to height:
###1) can't detect anything  < .1 m modeled height tall
###2) There are 9 'missed' plants known to exist. Each plant is taller than .1 m;
#### 8/9 modeled heights are < .1 m.
###3) So, the imperfect detection seems both stochastically related to real height
###(smaller plants likely to be missed? maybe not enough info?) and deterministically related to modeled height
###(if the height is modeled as < .1 m, the plant cannot be detected).
###This is a little trickier to think about...

###A global size distribution for the modeled height is left-truncated at .1. 
###I'm not sure why we would want to estimate this, or why we would want to 
###estimate the abundance of shrubs greater than some modeled height.

###A global size distribution for the 'true' height may not be left-truncated at .1 m:
###in theory, a modeled height could be > .1 m while the true height is < .1 m. 
###I think we do want to know N shrubs greater than 10 cm tall. But these may fall
###below the size threshold. Anyway, some edits below reflect this.


###Recall:
###y=0, 1 depending on on whether individual shrub is seen or not
###z is sampled; 1=shrub truly exists, 0 = shrub does not truly exist
###most z's are unknown initially (NA).
###Two exceptions: there will be some shrubs we know are falsely detected
###y=1, but z=0
###There will be some shrubs that exist that we know are not detected
###z=1, y=0. 
###Let's just assume using nimble to start.
library(nimble)
library(nimble)
Mod<-nimbleCode({

for (b in 1:nfx){
  beta[b]~dnorm(0, sd=1) ###describe intensity of point process
}

p11~dunif(0, 1)

###False positive parameters
p10~dunif(0, 1)


###size hyper-parameters. Reminder that these need to be checked!
for (c in 1:2){
shape[c]~dunif(0, 50)
rate[c]~dunif(0, 50)
}
###need to someway to relate true/modeled height
kappa0~dnorm(0, sd=1)
kappa1~dnorm(0, sd=1) ###maybe this needs to be informative
sigma~dgamma(1, 5)
#sigma<-dnorm(0, sd=1) T(0, ) truncation not working in v0.12?
###maybe want a different prior here for mixing

for (j in 1:nPix){
  mu[j]<-exp(inprod(beta[1:nfx], X[j, ]))  
  ##above implies a poisson point process. Potentially not reasonable.....
  probs[j]<-mu[j]/EN
}

EN<-sum(mu[1:nPix]) ###Expected total abundance
psi<-EN/M ###M is total number of known to exist, observed, and augemented shrubs

for (i in 1:M){
  z[i]~dbern(psi) ###does the shrub exist
  s[i]~dcat(probs[1:nPix]) ###where--which pixel--does the shrub exist
  ###if you trace the s[i] (where z[i]=1), you have estimates of the realized
  ###or actual, on ground (vs. expected) distribution of shrubs
  
  ###don't know the size of the unseen shrubs, so this has to be sampled...
  ###simple example, we just treat size as a gamma RV
  ###if s[i] describes the location of points
  ###size[i] is analagous to a 'mark' for the points

  
  ht_true[i]~dgamma(shape[z[i]+1], rate[z[i]+1]) #T(.1, ) 
  ###weird error in v0.12 with T. But not sure we want to truncate this? 
  #Do we think every possibly real (observed) shrub in the dataset is truly > .1 m tall?
  height[i]~dnorm(kappa0+kappa1*ht_true[i], sd=sigma)
  
  ###Slightly weird part. Can't be detected at all if modeled height
  ###is less than .1m.
  y[i]~dbern(step(height[i]-.1)*(z[i]*p11+(1-z[i])*p10))
}

N<-sum(z[1:M])
N2<-sum(z[1:2846]) ###how many of the observed plants exist...
})

###Couple places for expansion.
###Could consider any number of different discrete point processes:
###e.g., could consider terms or parmaeterizations for mu with more or less dispersion.
###Could consider spatial variation in size[i]...e.g.:
###size[i]~dgamma(shape[s[i]], rate[s[i]])
###where shape and rate vary across nPix in some fashion.
###Common to use a log link for gamma regression, something like...
###linear_predictor[j]<-inprod(delta[1:nfz], X[j, ])
###size[i]~dgamma(shape, shape / exp(linear_predictor[s[i]]))
###Above, mean--but not dispersion vary by pixel. Probably want the dispersion 
###to change too. 


###Anyway, on to the rest of the model parts needed to put this together.

Constants<-list(nfx=2, M=nrow(SpDat)+1000, nPix=nrow(CellDat)) ###M includes augmented individuals.
Data<-list(ht_true=c(SpDat$ht_true, rep(NA, 1000)), height=c(SpDat$height, rep(NA, 1000)), 
           X=cbind(rep(1, nrow(CellDat)), as.numeric(scale(CellDat$dist))),
           y=c(SpDat$y, rep(0, 1000)), z=c(SpDat$z, rep(NA, 1000)),
           s=c(SpDat$cellid, rep(NA, 1000)))
###Next initial values. We might have to be careful here...without
###thinking hard about it, might be possible for the sampler (or us)
###to set up initial values which really fuck things up.
###What could be initialized: beta, alpha, kappa, delta, z...
###height, ht_true...
###for beta, alpha, other params...can provide whatever is sensible.

###for something also provided as data (z, height, height true)
###if wanting to provide inits, need to do something like this...
zIn<-rep(NA, length(Data$z))
zIn[is.na(Data$z)]<-1
##above, we initialize every shrub as existing (not strictly neccessary). Note, we do not initialize
##the non-NA values passed as data! If we know they exist or not, don't want to pass an initial value.

Inits<- list(z=zIn)

Shrub<- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                        data=Data, inits = Inits)

ShrubConf<- configureMCMC(Shrub, monitors = c("N", "N2", 'beta', 'kappa0',
                                              'kappa1', 'p11', 'p10'))
##monitors--what things to trace


Rmcmc<-buildMCMC(ShrubConf)
compMCMC <- compileNimble(Rmcmc, Shrub)

samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=5000, nburnin=3000, thin=5, 
               nchains=1)

###Will obviously want to check/play around with chain length, add chains, etc. 
###This takes a few minutes on a dinky laptop.
           
