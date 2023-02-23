
###Note, probably a bug/error in here or two
SimDrone<-function(nsims, ncells, alphaTrue, betaTrue, alphaFalse, betaFalse, shape_false, rate_false,
                   shape_true, rate_true, delta, nVal){
N_True_Val<-rep(NA, nsims)
N_False_Val<-rep(NA, nsims)
N_True<-rep(NA, nsims)
N_est<-rep(NA, nsims)
beta_est<-matrix(NA, nsims, 2)
for (s in 1:nsims){
  print(s)
  ###create landscape
  grid<-expand.grid(1:sqrt(ncells[s]), 1:sqrt(ncells[s]))
  colnames(grid)<-c("X", "Y")
  grid$cell<-1:nrow(grid)
  BasisX<-splines::bs(grid$X, df = 5, degree = 3, intercept = FALSE)
  BasisY<-splines::bs(grid$Y, df = 5, degree = 3, intercept = FALSE)
  pBurn<-plogis(BasisX %*% rnorm(5, 0, 1.5)+BasisY %*% rnorm(5, 0, 1.5))
  pBurn<-pBurn-mean(pBurn)+.5
  grid$Burn<-rbinom(nrow(grid), 1, pBurn)
  grid$Burn[is.na(grid$Burn)]<-0
  
  ##allocate individuals
  grid$TrueDens<-exp(alphaTrue[s]+betaTrue[s]*grid$Burn)
  grid$TrueN<-rpois(nrow(grid), grid$TrueDens)
  ###Want to store this...
  N_True[s]<-sum(grid$TrueN)
  
  grid$FalseDens<-exp(alphaFalse[s]+betaFalse[s]*grid$Burn)
  grid$FalseN<-rpois(nrow(grid), grid$FalseDens)
  
  ###"distribute" the individuals--cast the N into an individual data frame
  FalseGuys<-data.frame(ID=1:sum(grid$FalseN), cellID=rep(grid$cell, grid$FalseN))
  TrueGuys<-data.frame(ID=1:sum(grid$TrueN), cellID=rep(grid$cell, grid$TrueN))
  
  ###generate heights
  FalseGuys$height<-rgamma(nrow(FalseGuys), shape_false[s], rate_false[s])
  FalseGuys<-FalseGuys[FalseGuys$height >= 0.1,]
  TrueGuys$height<-rgamma(nrow(TrueGuys), shape_true[s], rate_true[s])
  
  ###simulate detections
  FalseGuys$y<-1
  for (i in 1:nrow(TrueGuys)){
    TrueGuys$y[i]<-ifelse(TrueGuys$height[i]<.1, 0, rbinom(1, 1, plogis(delta[s, 1]+delta[s, 2]*TrueGuys$height[i])))
  }
  
  ###Validation process...
  Validation<-sample(grid$cell, nVal[s])
  
  N_True_Val[s]<-sum(TrueGuys$cellID %in% Validation)
  ###number of confirmed false guys
  N_False_Val[s]<-sum(FalseGuys$cellID %in% Validation)
  ###Number of burned cells validated
  sum(grid$Burn[Validation])
  
  `%notin%` <- Negate(`%in%`)
  NonVals<-rbind(FalseGuys[FalseGuys$cellID %notin% Validation, ], TrueGuys[TrueGuys$cellID %notin% Validation, ])
  Vals<-rbind(FalseGuys[FalseGuys$cellID %in% Validation, ], TrueGuys[TrueGuys$cellID %in% Validation, ])
  
  NonVals$z<-NA
  Vals$z<-c(rep(0, sum(FalseGuys$cellID %in% Validation)), rep(1, sum(TrueGuys$cellID %in% Validation)))
  
  ###We don't see any of the NonVals if modeled height is <.1 m
  NonVals<-NonVals[NonVals$height >=.1,]
  
  library(nimble)
  Mod<-nimbleCode({
    
    for (b in 1:nfx){
      for (c in 1:2){
        beta[b, c]~dnorm(0, sd=1) ###describe intensity of point processes
      }
    }
    
    alpha_0~dlogis(0, 1)
    alpha_1~dlogis(0, 1)
    
    logit(p10)~dlogis(0, 1) 
    for (c in 1:2){
      shape[c]~dunif(0, 10)
      rate[c]~dunif(0, 50)
    }
    
    ###nimble does not accept vectors of indices...
    ###So...must split the cells into sets
    
    for (j in 1:nNaive){
      mu[j, 2]<-exp(inprod(beta[1:nfx, 2], XNaive[j, 1:nfx]))  
      probs[j, 2]<-mu[j, 2]/EN[2, 1]
      mu[j, 1]<-exp(inprod(beta[1:nfx, 1], XNaive[j, 1:nfx]))  
      probs[j, 1]<-mu[j, 1]/EN[1, 1]
    }
    
    for (j in 1:nVal){
      mu2[j, 2]<-exp(inprod(beta[1:nfx, 2], XVal[j, 1:nfx]))  
      probs2[j, 2]<-mu2[j, 2]/EN[2, 2]
      mu2[j, 1]<-exp(inprod(beta[1:nfx, 1], XVal[j, 1:nfx]))  
      probs2[j, 1]<-mu2[j, 1]/EN[1, 2]
    }
    
    
    EN[1, 1]<-sum(mu[1:nNaive, 1]) ###Expected  abundance of fakes
    EN[2, 1]<-sum(mu[1:nNaive, 2]) ###Expected abundance of reals in unvalidated cells
    EN[1, 2]<-sum(mu2[1:nVal, 1]) ###Expected abundance of fakes in validated cells
    EN[2, 2]<-sum(mu2[1:nVal, 2]) ###expected abundance of reals in validated cells
    psi<-sum(EN[2, 1:2])/(M+V) ###M and V is total number of known to exist, observed, and augemented shrubs
    
    for (i in 1:M){ ###M this is the unverified data
      z[i]~dbern(psi) ###does the shrub exist
      cl[i]<-z[i]+1
      s[i]~dcat(probs[1:nNaive, cl[i]]) ###where--which pixel--does the shrub exist
      height[i]~dgamma(shape[cl[i]], rate[cl[i]])
      logit(p11[i])<-alpha_0+alpha_1*height[i]
      y[i]~dbern(step(height[i]-.1)*(z[i]*p11[i]+(1-z[i])*p10))
    }
    for (i in 1:V){ ###M this is the verified data
      z2[i]~dbern(psi) ###does the shrub exist
      class2[i]<-z2[i]+1
      s2[i]~dcat(probs2[1:nVal, class2[i]]) ###where--which pixel--does the shrub exist
      height2[i]~dgamma(shape[cl[i]], rate[cl[i]])
      logit(p12[i])<-alpha_0+alpha_1*height2[i]
      y2[i]~dbern(step(height2[i]-.1)*(z2[i]*p12[i]+(1-z2[i])*p10))
    }
    N<-sum(z[1:M])+sum(z2[1:V])
    ###how many of the observed plants exist...
  })
  
  NonValidatedCells<-grid$cell[grid$cell %notin% Validation]
  NonVals$s<-NA
  for (i in 1:nrow(NonVals)){
    NonVals$s[i]<-which(NonValidatedCells==NonVals$cellID[i])
  }
  
  Vals$s<-NA
  for (i in 1:nrow(Vals)){
    Vals$s[i]<-which(Validation==Vals$cellID[i])
  }
  
  
  Constants<-list(M=nrow(NonVals)+1000, V=nrow(Vals), 
                  nfx=2, nNaive=nrow(grid)-length(Validation), nVal=length(Validation),
                  XNaive=cbind(rep(1, nrow(grid)-length(Validation)), grid$Burn[-Validation]),
                  XVal=cbind(rep(1, length(Validation)), grid$Burn[Validation]), 
                  class2=Vals$z+1)
  
  Data<-list(y=c(NonVals$y, rep(0, 1000)), s=c(NonVals$s, rep(NA, 1000)),
             z2=Vals$z, y2=Vals$y, s2=Vals$s, height=c(NonVals$height, rep(NA, 1000)), 
  hinit<-c(rep(NA, nrow(NonVals)), rgamma(1000, 2, 10)) 
  Zinit<-c(rep(1, nrow(NonVals)), rep(0, 1000))
  Cinit<-Zinit+1
  Inits <- list(z=Zinit, cl=Cinit, height=hinit, psi=0.5, p10=0.75, alpha_0=0, alpha_1=0, beta=matrix(0, 2, 2), shape=rep(5, 2), rate=rep(10, 2))
             
  Shrub <- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                                  data=Data, inits = Inits)
             
  mcmcConf <- configureMCMC(Shrub, monitors = c("alpha_0", "alpha_1",  "N", "beta", 
          "psi", "p10", "shape", "rate"), useConjugacy=FALSE)
             
  Rmcmc<-buildMCMC(mcmcConf)
             
             
  compMCMC <- compileNimble(Rmcmc, Shrub)
             
  samps<-runMCMC(mcmc = compMCMC$Rmcmc,
                            niter=8000, nburnin=5000, thin=1, 
                            nchains=1)
             
  
  ##extract info from model...
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  N_est[s]<-mode(samps[, "N"]
  ###probably should use 90% CI unless running for a long time....
  N_lci[s]<-quantile(samps[, "N"], probs=.05)
  N_uci[s]<-quantile(samps[, "N"], probs=.95)
  beta_est[s, ]<-c(mean(samps[, "beta[1, 2]"]), mean(samps[, "beta[2, 2]"]))
  
  ###anything else we want to store?
} ###end sims loop
return(N_est=N_est, N_lci=N_lci, N_uci=N_uci, beta_est=beta_est, N_True_Val=N_True_Val,
       N_False_Val=N_False_Val, N_True=N_True)
}

}