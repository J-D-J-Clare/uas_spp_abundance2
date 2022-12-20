

###create grid, covariates, individuals
grid<-expand.grid(1:40, 1:40)
colnames(grid)<-c("X", "Y")
grid$cell<-1:nrow(grid)
BasisX<-splines::bs(grid$X, df = 5, degree = 3, intercept = FALSE)
BasisY<-splines::bs(grid$Y, df = 5, degree = 3, intercept = FALSE)
pBurn<-plogis(BasisX %*% rnorm(5, 0, 1.5)+BasisY %*% rnorm(5, 0, 1.5))
pBurn<-pBurn-mean(pBurn)+.5
grid$Burn<-rbinom(nrow(grid), 1, pBurn)
grid$Burn[is.na(grid$Burn)]<-0
#ggplot(grid, aes(X, Y, Burn))+geom_point(aes(color=factor(Burn)))+scale_color_viridis_d()

grid$TrueDens<-exp(1-1.5*grid$Burn)
#ggplot(grid, aes(X, Y, TrueDens))+geom_point(aes(color=TrueDens))+scale_color_viridis_b()
grid$TrueN<-rpois(nrow(grid), grid$TrueDens)
#ggplot(grid, aes(X, Y, TrueN))+geom_point(aes(color=TrueN))+scale_color_viridis_b()


grid$FalseDens<-exp(-.5+.5*grid$Burn)
#ggplot(grid, aes(X, Y, FalseDens))+geom_point(aes(color=FalseDens))+scale_color_viridis_b()
grid$FalseN<-rpois(nrow(grid), grid$FalseDens)
library(ggplot2)
#ggplot(grid, aes(X, Y, FalseN))+geom_point(aes(color=FalseN))+scale_color_viridis_b()
a<-ggplot(grid, aes(X, Y, Burn))+geom_point(aes(color=factor(Burn)))+scale_color_viridis_d()
b<-ggplot(grid, aes(X, Y, TrueDens))+geom_point(aes(color=TrueDens))+scale_color_viridis_b()
c<-ggplot(grid, aes(X, Y, TrueN))+geom_point(aes(color=TrueN))+scale_color_viridis_b()
e<-ggplot(grid, aes(X, Y, FalseDens))+geom_point(aes(color=FalseDens))+scale_color_viridis_b()
f<-ggplot(grid, aes(X, Y, FalseN))+geom_point(aes(color=FalseN))+scale_color_viridis_b()

gridExtra::grid.arrange(a, b, c, e, f, ncol=3)





###"distribute" the individuals--cast the N into an individual data frame
FalseGuys<-data.frame(ID=1:sum(grid$FalseN), cellID=rep(grid$cell, grid$FalseN))
TrueGuys<-data.frame(ID=1:sum(grid$TrueN), cellID=rep(grid$cell, grid$TrueN))

###"False" plants have only a modeled size
### This could be a reason height has not been informative...
### false guys don't really follow a gamma distribution, but a weird truncated distribution.
### would need to think about how to model this.
FalseGuys$realHt<-NA
FalseGuys$height<-rgamma(nrow(FalseGuys), 1, 25)+.1

###True plants have a real height, and then a modeled height.
###presumably, the true height inversely relates to the density within the cell?
###or, it could be that greater density increases the dispersion--more seedlings with more mature plants?
grid$ESize<-exp(.25-.75*grid$TrueN-.5*grid$Burn)
hist(grid$ESize) ###note, far right reflects N=0, so these plants do not exist
grid$shape<-simstudy::gammaGetShapeRate(grid$ESize, .25)$shape
grid$rate<-simstudy::gammaGetShapeRate(grid$ESize, .25)$rate
TrueGuys$realHt<-rgamma(nrow(TrueGuys), shape=grid$shape[TrueGuys$cellID], rate=grid$rate[TrueGuys$cellID])
hist(TrueGuys$realHt)

###the height is modeled with some error...
TrueGuys$height<-rnorm(nrow(TrueGuys), -.05+.9*TrueGuys$realHt, .07)
plot(TrueGuys$realHt, TrueGuys$height)
sum(TrueGuys$realHt>.1) ###I guess this is the true abundance for this problem?

###Detection process. The false guys are all 'seen'.
###Some of the true guys will be missed:
###deterministically if their modeled height is less than whatever threshold;
###stochastically depending upon size. 
FalseGuys$y<-1
for (i in 1:nrow(TrueGuys)){
  TrueGuys$y[i]<-ifelse(TrueGuys$height[i]<.1, 0, rbinom(1, 1, plogis(-3+15*TrueGuys$realHt[i])))
}
###Just to give a sense of the association
plot(TrueGuys$realHt, plogis(-3+15*TrueGuys$realHt)) 
plot(TrueGuys$realHt, TrueGuys$y)
plot(TrueGuys$height, TrueGuys$y)

###Could also add some simulated measure of confidence in the plant
###following what might be produced by a CNN.




###Haven't created this yet, but for all true guys, the generating z[i]=1, for all false guys, z[i]=0

###Some number of 'grid cells' are validated in some way in the field.
###Many ways to consider a sampling scheme for the validation.

###Anyway, here, let's imagine 50 cells randomly sampled.
Validation<-sample(grid$cell, 50)

library(raster)
plot(rasterFromXYZ(grid[, c(1, 2, 4)]))
points(grid[Validation, 1:2])

###number of confirmed true guys
sum(TrueGuys$cellID %in% Validation)
###number of confirmed false guys
sum(FalseGuys$cellID %in% Validation)
###Number of burned cells validated
sum(grid$Burn[Validation])

`%notin%` <- Negate(`%in%`)
NonVals<-rbind(FalseGuys[FalseGuys$cellID %notin% Validation, ], TrueGuys[TrueGuys$cellID %notin% Validation, ])
Vals<-rbind(FalseGuys[FalseGuys$cellID %in% Validation, ], TrueGuys[TrueGuys$cellID %in% Validation, ])

NonVals$z<-NA
Vals$z<-ifelse(is.na(Vals$realHt), 0, 1)
#Vals$y[is.na(Vals$realHt)]<-1

###We don't see any of the NonVals if modeled height is <.1 m
NonVals<-NonVals[NonVals$height >=.1,]

###We can't "validate" the Vals if the True height < .1m
Vals<-Vals[-which(Vals$z==1 & Vals$realHt <=.1),]


###define the model
###For simplicity, may want to ignore size effects?

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
  kappa0~dnorm(0, sd=1)
  kappa1~dnorm(0, sd=1) ###maybe this needs to be informative
  sigma~dgamma(1, 5)
  
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
    ht_true[i]~T(dgamma(shape[cl[i]], rate[cl[i]]), .1, )
    #ht_true[i]~T(dgamma(shape[z[i]+1], rate[z[i]+1]), .1, ) 
    #Not sure we want to truncate this? 
    height[i]~dnorm(kappa0+kappa1*ht_true[i], sd=sigma)
    logit(p11[i])<-alpha_0+alpha_1*ht_true[i]
    y[i]~dbern(step(height[i]-.1)*(z[i]*p11[i]+(1-z[i])*p10))
  }
  for (i in 1:V){ ###M this is the verified data
    z2[i]~dbern(psi) ###does the shrub exist
    #class2[i]<-z2[i]+1
    s2[i]~dcat(probs2[1:nVal, class2[i]]) ###where--which pixel--does the shrub exist
    ht_true2[i]~T(dgamma(shape[class2[i]], rate[class2[i]]), .1, )
    ##ht_true[i]~T(dgamma(shape[z[i]+1], rate[z[i]+1]), .1, ) 
    height2[i]~dnorm(kappa0+kappa1*ht_true2[i], sd=sigma)
    logit(p12[i])<-alpha_0+alpha_1*ht_true2[i]
    y2[i]~dbern(step(height[i]-.1)*(z2[i]*p12[i]+(1-z2[i])*p10))
  }
  N<-sum(z[1:M])+sum(z2[1:V])
  ###how many of the observed plants exist...
})

###Right, need to re-index the cells
NonValidatedCells<-grid$cell[grid$cell %notin% Validation]
NonVals$s<-NA
for (i in 1:nrow(NonVals)){
  NonVals$s[i]<-which(NonValidatedCells==NonVals$cellID[i])
}

Vals$s<-NA
for (i in 1:nrow(Vals)){
  Vals$s[i]<-which(Validation==Vals$cellID[i])
}


Constants<-list(M=nrow(NonVals+1000), V=nrow(Vals), 
                nfx=2, nNaive=nrow(grid)-length(Validation), nVal=length(Validation),
                XNaive=cbind(rep(1, nrow(grid)-length(Validation)), grid$Burn[-Validation]),
                XVal=cbind(rep(1, length(Validation)), grid$Burn[Validation]), 
                class2=Vals$z+1)

Data<-list(y=c(NonVals$y, rep(0, 1000)), s=c(NonVals$s, rep(NA, 1000)),
           z2=Vals$z, y2=Vals$y, s2=Vals$s, height=c(NonVals$height, rep(NA, 1000)), 
           ht_true=rep(NA, nrow(NonVals)+1000), height2=Vals$height,
           ht_true2=Vals$realHt)
#validated=Validation, nonvalidated=grid$cell %notin% Validation) 

Zinit<-c(rep(1, nrow(NonVals)), rep(0, 1000))
Cinit<-Zinit+1
Inits <- list(z=Zinit, cl=Cinit, beta=matrix(0, 2, 2), ht_true=runif(nrow(NonVals)+1000, .1, .3))

Shrub <- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                     data=Data, inits = Inits)

mcmcConf <- configureMCMC(Shrub, monitors = c("alpha_0", "alpha_1", "p10", "N", "beta", 
                                              "kappa0", "kappa1", "sigma"))

Rmcmc<-buildMCMC(mcmcConf)


compMCMC <- compileNimble(Rmcmc, Shrub)

samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=5000, nburnin=3000, thin=1, 
               nchains=1)

###needs to run longer...some parameters in the right ballpark, at least. 
###I kind of suspect that sampling the latent existence states
###and having each state have it's own "distribution" parameters
###is tricky. N not looking so good.

###discussion for meeting--

###1) Primarily, would be good to assess data generation process.
### Anything missing or non-desirable?

###2) This is not yet wrapped up in a simulation function--easy to do--
###where validation sample, strategy, other parameters could be tweaked.
###Be good to start to hit on a range of variation for the parameters of interest.

###3) Kind of a side note, but may want to balance complexity and speed. 
### Size seems like a big deal in practice, but adds a lot of moving parts.
### Can do this much more quickly as an N-mixture sort of model, but
### can't have individual variation in p. 