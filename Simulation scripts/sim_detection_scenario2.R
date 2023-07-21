nsims=200
nnonval=1400
alphaTrue=1.5 
betaTrue=.7
alphaFalse=0.25 
betaFalse1=1
betaFalse2=-.5 
p1=.7 ###ground-based p
p2=.5 ###air-based p
nVal=200 

N_True<-rep(NA, nsims) ###sum of N over each "unit"

XNaive<-array(NA, dim=c(nnonval, 3, nsims))
XVal<-array(NA, dim=c(nVal, 3, nsims))
XNaive[, 1, ]<-1
XVal[, 1, ]<-1

###data includes:
###c_air (aerial counts that are mixes of true positives and false positives in Naive cells)
cstar_air<-matrix(NA, nnonval, nsims)


###In validated cells, there are:
##a) two independent ground counts
##b) aerial counts that intersect with the ground counts (TP's)
##c) aerial counts that do not intersect with the ground counts (TPs, FP's)
##hence, the count is multivariate, and perhaps expressed as a 
##multinomial....



#cTP_air<-matrix(NA, nVal, nsims)
#cstar_air<-matrix(NA, nVal, nsims)
#c_ground<-array(NA, dim=c(nVal, nsims, 2)) ###assuming 2 replicate counts
n_ground<-matrix(NA, nVal, nsims) ###holding this for now--this is the total count on the ground
y<-array(NA, dim=c(nVal, nsims, 7))

###the N[j] individuals will be split into 8 boxes...
###1 detected first ground survey, not second, not in air
###2 not detected first ground survey, detected second, not in air
###3 detected first ground survey, detected second, not in air
###4 detected first ground survey, not-detected second, detected in air
###5 not detected first ground survey, detected second, detected in air
###6 detected first ground survey, detected second, detected in air
###7 not detected on ground, detected in air
###8 not observed at all
probvec<-rep(NA, 8)
probvec[1]<-p1*(1-p1)*(1-p2)
probvec[2]<-(1-p1)*p1*(1-p2)
probvec[3]<-p1*p1*(1-p2)
probvec[4]<-p1*(1-p1)*p2
probvec[5]<-(1-p1)*p1*p2
probvec[6]<-p1*p1*p2
probvec[7]<-(1-p1)*(1-p1)*p2
probvec[8]<-(1-p1)*(1-p1)*(1-p2)

###n_ground equals sum of n's in 1:6<-
###c_ground[,,1] equals sum of 1, 3, 4, 6
###c_ground[,,2] equals sum of 2, 3, 5, 6
###c_tp_air equals sums 4:6
###c_star air equals fp's plus 7



for (i in 1:nsims){
  XNaive[, 2, i]<-rnorm(nnonval)
  XNaive[, 3, i]<-XNaive[, 2, i]^2
  XVal[, 2, i]<-rnorm(nVal)
  XVal[, 3, i]<-XVal[, 2, i]^2
  
  TrueMuNaive<-exp(alphaTrue+betaTrue*XNaive[, 2, i])
  TrueNNaive<-rpois(nnonval, TrueMuNaive)
  
  TrueMuVal<-exp(alphaTrue+betaTrue*XVal[, 2, i])
  TrueNVal<-rpois(nVal, TrueMuVal)
  
  ###Want to store this...
  N_True[i]<-sum(TrueNNaive, TrueNVal)
  
  ###raw aerial counts where there is no ground sampling
  cstar_air[, i]<-rpois(nnonval, exp(alphaFalse+betaFalse1*XNaive[, 2, i]+betaFalse2*XNaive[, 3, i]))+
    rbinom(nnonval, TrueNNaive, p2)
  
  ###now--a little more involved--work through the cells with ground counts also
  ### here's the false aerial counts in the cells with ground sampling
  c_air_fp<-rpois(nVal, exp(alphaFalse+betaFalse1*XVal[, 2, i]+betaFalse2*XVal[, 3, i]))
  for (x in 1:nVal){
    tmp<-as.vector(rmultinom(1, TrueNVal[x], probvec))
    ###n_ground equals sum of n's in 1:6
    ###c_star air equals fp's plus 7
    y[x, i, 1:7]<-tmp[1:7]
    n_ground[x, i]<-sum(tmp[1:6])
    y[x, i, 7]<-y[x, i, 7]+c_air_fp[x]
  }
}




library(nimble)
dCompoundNMix<-nimbleFunction(
  run = function(x = integer(0), N = integer(0), p1=double(0), lamFP=double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    ll1<-nimRep(0, x+1)
    ll2<-nimRep(0, x+1)
    for (i in 1:(x+1)) {
      ll1[i]<-dbinom(i-1, N, p1)
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
    
    p1[s]~dunif(0, 1) ###pr shrub detected via single ground count
    p2[s]~dunif(0, 1) ###pr shrub detected from air
    p_g_n[s]<-(p1[s]^2)+2*(p1[s]*(1-p1[s])) ##pr shrub detected on at least 1 ground count
    
    pi[s, 1]<-(p1[s]*(1-p1[s])*(1-p2[s]))/p_g_n[s] # NTP[j, s]
    pi[s, 2]<-((1-p1[s])*p1[s]*(1-p2[s]))/p_g_n[s]
    pi[s, 3]<-(p1[s]*p1[s]*(1-p2[s]))/p_g_n[s]
    pi[s, 4]<-(p1[s]*(1-p1[s])*p2[s])/p_g_n[s]
    pi[s, 5]<-((1-p1[s])*p1[s]*p2[s])/p_g_n[s]
    pi[s, 6]<-(p1[s]*p1[s]*p2[s])/p_g_n[s]
    
    
    
    for (j in 1:nNaive){
      mu[j, 2, s]<-exp(inprod(beta[1:nfx[2], 2, s], XNaive[j, 1:nfx[2], s]))
      mu[j, 1, s]<-exp(inprod(beta[1:nfx[1], 1, s], XNaive[j, 1:nfx[1], s]))  ##fp
      N[j, s]~dpois(mu[j, 2, s])
      cstar_air[j, s]~dCompoundNMix(N=N[j,s], p=p2[s], lamFP=mu[j, 1, s])
    }
    
    for (j in 1:nVal){ ##counts in all validated cells
      mu2[j, 2, s]<-exp(inprod(beta[1:nfx[2], 2, s], XVal[j, 1:nfx[2], s]))
      mu2[j, 1, s]<-exp(inprod(beta[1:nfx[1], 1, s], XVal[j, 1:nfx[1], s])) ##fp
      NTP[j, s]~dpois(mu2[j, 2, s])
      n_ground[j, s]~dbin(p_g_n[s], NTP[j, s])
      y[j, s, 1:6]~dmulti(size=n_ground[j, s], prob=pi[s, 1:6])
      y[j, s, 7]~dCompoundNMix(N=NTP[j,s], p=(1-p1[s])^2*p2[s], lamFP=mu2[j, 1, s])
    }
    
    
    
    ###guess at the number of shrubs that exist
    EN[s]<-sum(mu[1:nNaive, 2, s])+sum(mu2[1:nVal, 2, s])
    TotalN[s]<-sum(N[1:nNaive, s])+sum(NTP[1:nVal, s])
  }
})


Constants<-list(nsims=nsims, nfx=c(3, 2), nNaive=nnonval, nVal=nVal)

Data<-list(y=y, n_ground=n_ground, cstar_air=cstar_air,
           XNaive=XNaive,
           XVal=XVal)

Inits <- list(p1=rep(.5, nsims), p2=rep(.5, nsims), beta=array(0, dim=c(3, 2, nsims)),
              NTP=n_ground, N=cstar_air)

Shrub <- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                     data=Data, inits = Inits)

mcmcConf <- configureMCMC(Shrub, monitors = c("p1", "p2", "beta", 
                                              "EN",  "TotalN"), useConjugacy=FALSE)

Rmcmc<-buildMCMC(mcmcConf)


compMCMC <- compileNimble(Rmcmc, Shrub)

samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=10000, nburnin=5000, thin=5, 
               nchains=3)
