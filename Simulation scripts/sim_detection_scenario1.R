# Simulations under: nnonval and nVal as {1400/200, 1450/150, 1500/100, 1550/50}

nsims=200
nnonval=1400
alphaTrue=1.5 
alphaFalse=0.25 
betaFalse1=1
betaFalse2=-.5 
p=.5 ###air-based p
nVal=200 ###cells

N_True<-rep(NA, nsims) ###sum of N over each "unit"

XNaive<-array(NA, dim=c(nnonval, 3, nsims))
XVal<-array(NA, dim=c(nVal, 3, nsims))
XNaive[, 1, ]<-1
XVal[, 1, ]<-1

###data here includes:
###NTP (ground-verified abundance)
###cTP_air (number of aerial detections "known" to be TP in verified cells)
###cFP_air (number of aerial detections "known" to be FP in verified cells)
###cstar_air (aerial detections that could be TP, could be FP at sites lacking ground sampling)


cstar_air<-matrix(NA, nnonval, nsims)
NTP<-matrix(NA, nVal, nsims)
cTP_air<-matrix(NA, nVal, nsims)
cFP_air<-matrix(NA, nVal, nsims)


for (i in 1:nsims){
  XNaive[, 2, i]<-rnorm(nnonval)
  XNaive[, 3, i]<-XNaive[, 2, i]^2
  XVal[, 2, i]<-rnorm(nVal)
  XVal[, 3, i]<-XVal[, 2, i]^2
  
  TrueMuNaive<-exp(alphaTrue+betaTrue*XNaive[, 2, i])
  TrueNNaive<-rpois(nnonval, TrueMuNaive)
  
  TrueMuVal<-exp(alphaTrue+betaTrue*XVal[, 2, i])
  NTP[, i]<-rpois(nVal, TrueMuVal)
  
  ###Want to store this...
  N_True[i]<-sum(TrueNNaive, NTP[, i])
  
  ###raw aerial counts where there is no ground sampling
  cstar_air[, i]<-rpois(nnonval, exp(alphaFalse+betaFalse1*XNaive[, 2, i]+betaFalse2*XNaive[, 3, i]))+
    rbinom(nnonval, TrueNNaive, p)
  
  ###data in verified cells
  cFP_air[, i]<-rpois(nVal, exp(alphaFalse+betaFalse1*XVal[, 2, i]+betaFalse2*XVal[, 3, i]))
  cTP_air[, i]<-rbinom(nVal, NTP[, i], p)
}




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
    else return(exp(logProb)) })



Mod<-nimbleCode({
  for (s in 1:nsims){
    for (c in 1:2){
      for (b in 1:nfx[c]){
        beta[b, c, s]~dnorm(0, sd=1) ###describe intensity of point processes
      }
    }
    
    p[s]~dunif(0, 1) ###pr real shrub detected from air
    for (j in 1:nNaive){
      mu[j, 2, s]<-exp(inprod(beta[1:nfx[2], 2, s], XNaive[j, 1:nfx[2], s]))
      mu[j, 1, s]<-exp(inprod(beta[1:nfx[1], 1, s], XNaive[j, 1:nfx[1], s]))  
      N[j, s]~dpois(mu[j, 2, s])
      cstar_air[j, s]~dCompoundNMix(N=N[j,s], p=p[s], lamFP=mu[j, 1, s])
    }
    
    for (j in 1:nVal){ ##counts in all validated cells
      mu2[j, 2, s]<-exp(inprod(beta[1:nfx[2], 2, s], XVal[j, 1:nfx[2], s]))
      mu2[j, 1, s]<-exp(inprod(beta[1:nfx[1], 1, s], XVal[j, 1:nfx[1], s]))
      NTP[j, s]~dpois(mu2[j, 2, s])
      cTP_air[j, s]~dbin(p[s], NTP[j, s]) ###nTP is now the number of tp detects on ground
      cFP_air[j, s]~dpois(mu2[j, 1, s])
    }
    
    
    ###guess at the number of shrubs that exist
    EN[s]<-sum(mu[1:nNaive, 2, s])+sum(mu2[1:nVal, 2, s])
    TotalN[s]<-sum(N[1:nNaive, s])+sum(NTP[1:nVal, s])
  }
})


Constants<-list(nsims=nsims, nfx=c(3, 2), nNaive=nnonval, nVal=nVal)

Data<-list(cFP_air=cFP_air, cTP_air=cTP_air, NTP=NTP, cstar_air=cstar_air,
           XNaive=XNaive,
           XVal=XVal)

Inits <- list(p=rep(.5, nsims), beta=array(0, dim=c(3, 2, nsims)),
              N=cstar_air+3)

Shrub <- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                     data=Data, inits = Inits)

mcmcConf <- configureMCMC(Shrub, monitors = c("p", "beta", 
                                              "EN",  "TotalN"), useConjugacy=FALSE)

Rmcmc<-buildMCMC(mcmcConf)


compMCMC <- compileNimble(Rmcmc, Shrub)

samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=10000, nburnin=5000, thin=5, 
               nchains=3)