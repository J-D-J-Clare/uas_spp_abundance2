nsims=200
ncells=1600
alphaTrue=-.5 ###intercept--expecting about .5 real plants per "average" cell--can modify!
betaTrue=.7
alphaFalse=0 ###intercept--expecting about 1 false plant per average cell--can modify!
betaFalse1=1
betaFalse2=-.5 
p=.9
nVal=50 ###cells

N_True<-rep(NA, nsims)
N_True_Val<-rep(NA, nsims)
###number of confirmed false guys
N_False_Val<-rep(NA, nsims)
#EN<-rep(NA, nsims)
#betahat<-matrix(NA, 6, nsims)
#phat<-rep(NA, nsims)
nseen<-matrix(NA, ncells-nVal, nsims)
nFPval<-matrix(NA, nVal, nsims)
nTPval<-matrix(NA, nVal, nsims)


XNaive<-array(NA, dim=c(ncells-nVal, 3, nsims))
XVal<-array(NA, dim=c(nVal, 3, nsims))
XNaive[, 1, ]<-1
XVal[, 1, ]<-1

for (s in 1:nsims){
  grid<-expand.grid(1:sqrt(ncells), 1:sqrt(ncells))
  colnames(grid)<-c("X", "Y")
  grid$cell<-1:nrow(grid)
  BasisX<-splines::bs(grid$X, df = 5, degree = 3, intercept = FALSE)
  BasisY<-splines::bs(grid$Y, df = 5, degree = 3, intercept = FALSE)
  grid$Cov<-BasisX %*% rnorm(5, 0, 1.5)+BasisY %*% rnorm(5, 0, 1.5)
  grid$Cov<-grid$Cov-mean(grid$Cov)
  
  ##allocate individuals
  grid$TrueDens<-exp(alphaTrue+betaTrue*grid$Cov)
  grid$TrueN<-rpois(nrow(grid), grid$TrueDens)
  ###Want to store this...
  N_True[s]<-sum(grid$TrueN)
  
  grid$FalseDens<-exp(alphaFalse+betaFalse1*grid$Cov+betaFalse2*grid$Cov^2)
  grid$FalseN<-rpois(nrow(grid), grid$FalseDens)
  N_False<-sum(grid$FalseN)
  ###"distribute" the individuals--cast the N into an individual data frame
  FalseGuys<-data.frame(ID=1:sum(grid$FalseN), cellID=rep(grid$cell, grid$FalseN))
  TrueGuys<-data.frame(ID=1:sum(grid$TrueN), cellID=rep(grid$cell, grid$TrueN))
  
  ###simulate detections
  FalseGuys$y<-1
  for (i in 1:nrow(TrueGuys)){
    TrueGuys$y[i]<-rbinom(1, 1, p)
  }
  
  ###Validation process...
  Validation<-sample(grid$cell, nVal)
  
  XNaive[, 2:3, s]<-cbind(grid$Cov[-Validation], grid$Cov[-Validation]^2)
  XVal[, 2:3, s]<-cbind(grid$Cov[Validation], grid$Cov[Validation]^2)
  
  N_True_Val[s]<-sum(TrueGuys$cellID %in% Validation)
  ###number of confirmed false guys
  N_False_Val[s]<-sum(FalseGuys$cellID %in% Validation)
  ###Number of burned cells validated
  
  `%notin%` <- Negate(`%in%`)
  NonVals<-rbind(FalseGuys[FalseGuys$cellID %notin% Validation, ], TrueGuys[TrueGuys$cellID %notin% Validation, ])
  Vals<-rbind(FalseGuys[FalseGuys$cellID %in% Validation, ], TrueGuys[TrueGuys$cellID %in% Validation, ])
  Vals$z<-c(rep(0, sum(FalseGuys$cellID %in% Validation)), rep(1, sum(TrueGuys$cellID %in% Validation)))
  
  NonValidatedCells<-grid$cell[grid$cell %notin% Validation]
  NonVals$s<-NA
  for (i in 1:nrow(NonVals)){
    NonVals$s[i]<-which(NonValidatedCells==NonVals$cellID[i])
  }
  
  Vals$s<-NA
  for (i in 1:nrow(Vals)){
    Vals$s[i]<-which(Validation==Vals$cellID[i])
  }
  
  
  ###number of things seen in non validated cells
  n<-rep(0, nrow(grid)-length(Validation))
  for (i in 1:length(n)){
    n[i]<-length(which(NonVals$s==i))
  }
  
  nseen[,s]<-n
  
  nFP<-rep(0, length(Validation)) ###number of false things per validated cell
  nTP<-rep(0, length(Validation))  ###number of true things
  for (i in 1:length(nFP)){
    nFP[i]<-length(which(Vals$s==i & Vals$z==0))
    nTP[i]<-length(which(Vals$s==i & Vals$z==1))
    
  }
  nFPval[,s]<-nFP
  nTPval[,s]<-nTP
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
    
    logit(p[s])~dlogis(0, 1) ###pr real shrub detected from air
    for (j in 1:nNaive){
      mu[j, 2, s]<-exp(inprod(beta[1:nfx[2], 2, s], XNaive[j, 1:nfx[2], s]))
      mu[j, 1, s]<-exp(inprod(beta[1:nfx[1], 1, s], XNaive[j, 1:nfx[1], s]))  
      #n[j,s]~dpois(p[s]*mu[j, 2, s]+mu[j, 1, s]) ###nice thing about poisson--sums and thinning still marginally poisson
      NTP[j, s]~dpois(mu[j, 2, s])
      n[j, s]~dCompoundNMix(N=NTP[j,s], p=p[s], lamFP=mu[j, 1, s])
    }
    
    for (j in 1:nVal){ ##counts in all validated cells
      mu2[j, 2, s]<-exp(inprod(beta[1:nfx[2], 2, s], XVal[j, 1:nfx[2], s]))
      nTP[j, s]~dpois(mu2[j, 2, s])
      mu2[j, 1, s]<-exp(inprod(beta[1:nfx[1], 1, s], XVal[j, 1:nfx[1], s]))
      nFP[j, s]~dpois(mu2[j, 1, s])
    }
    
    ###validated plants...
    ###guess at the number of unvalided shrubs that exist
    EN[s]<-sum(mu[1:nNaive, 2, s])+sum(mu2[1:nVal, 2, s])
    N[s]<-sum(NTP[1:nNaive, s])+sum(nTP[1:nVal, s])
  }
})


Constants<-list(nsims=nsims, nfx=c(3, 2), nNaive=nrow(grid)-length(Validation), nVal=length(Validation))

Data<-list(n=nseen, nFP=nFPval, nTP=nTPval, 
           XNaive=XNaive,
           XVal=XVal)

Inits <- list(p=rep(0, nsims), beta=array(0, dim=c(3, 2, nsims)),
              NTP=matrix(1, Constants$nNaive, nsims))

Shrub <- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                     data=Data, inits = Inits)

mcmcConf <- configureMCMC(Shrub, monitors = c("p", "beta", 
                                              "EN",  "N"), useConjugacy=FALSE)

Rmcmc<-buildMCMC(mcmcConf)


compMCMC <- compileNimble(Rmcmc, Shrub)

samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=10000, nburnin=5000, thin=5, 
               nchains=1)

write.csv(samps, file="SampsVal50.csv")


params<-data.frame(cbind(N_True, N_True_Val, N_False_Val))
write.csv(params, file="ParamsVal50.csv")


