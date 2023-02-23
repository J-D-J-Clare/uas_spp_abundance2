

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

#Commenting out plotting commands
library(ggplot2)
ggplot(grid, aes(X, Y, Burn))+geom_point(aes(color=factor(Burn)))+scale_color_viridis_d()

grid$TrueDens<-exp(1-1.5*grid$Burn)
#ggplot(grid, aes(X, Y, TrueDens))+geom_point(aes(color=TrueDens))+scale_color_viridis_b()
grid$TrueN<-rpois(nrow(grid), grid$TrueDens)
#ggplot(grid, aes(X, Y, TrueN))+geom_point(aes(color=TrueN))+scale_color_viridis_b()


grid$FalseDens<-exp(.5+1*grid$Burn)
#ggplot(grid, aes(X, Y, FalseDens))+geom_point(aes(color=FalseDens))+scale_color_viridis_b()
grid$FalseN<-rpois(nrow(grid), grid$FalseDens)


ggplot(grid, aes(X, Y, FalseN))+geom_point(aes(color=FalseN))+scale_color_viridis_b()
a<-ggplot(grid, aes(X, Y, Burn))+geom_point(aes(color=factor(Burn)))+scale_color_viridis_d()
b<-ggplot(grid, aes(X, Y, TrueDens))+geom_point(aes(color=TrueDens))+scale_color_viridis_b()
c<-ggplot(grid, aes(X, Y, TrueN))+geom_point(aes(color=TrueN))+scale_color_viridis_b()
e<-ggplot(grid, aes(X, Y, FalseDens))+geom_point(aes(color=FalseDens))+scale_color_viridis_b()
f<-ggplot(grid, aes(X, Y, FalseN))+geom_point(aes(color=FalseN))+scale_color_viridis_b()
gridExtra::grid.arrange(a, b, c, e, f, ncol=3)





###"distribute" the individuals--cast the N into an individual data frame
FalseGuys<-data.frame(ID=1:sum(grid$FalseN), cellID=rep(grid$cell, grid$FalseN))
TrueGuys<-data.frame(ID=1:sum(grid$TrueN), cellID=rep(grid$cell, grid$TrueN))



### In reality, there is a modeled height and field-measured height...
### These are associated, but this association poses challenges.

### So, here's the deal instead:
### just give all true and false plants a modeled height.
### Those with a modeled height < .1 m can't be detected.

### For the false guys, this seems like a shifted gamma
### (Don't really care about "undetected" false phenomena)
###For simplicity, much easier to use a gamma.
##Would be good to know what the appropriate parms are here

FalseGuys$height<-rgamma(nrow(FalseGuys), 2, 20)

###Not gonna see any false guys shorter than 0.1m
FalseGuys<-FalseGuys[FalseGuys$height >= 0.1,]

###For true plants, this varies across space.
###below, smaller plants where there are more plants &
###smaller plants in burned areas.

###Not even sure we need to worry about this?
###Key question is whether to consider this in the model (post-simulation) or not?

grid$ESize<-exp(.5-.5*grid$TrueN-.5*grid$Burn)
hist(grid$ESize) ###note, far right reflects N=0, so these plants do not exist
grid$shape<-simstudy::gammaGetShapeRate(grid$ESize, .25)$shape
grid$rate<-simstudy::gammaGetShapeRate(grid$ESize, .25)$rate
TrueGuys$height<-rgamma(nrow(TrueGuys), shape=grid$shape[TrueGuys$cellID], rate=grid$rate[TrueGuys$cellID])

###Detection process. The false guys are all 'seen'.
###Some of the true guys will be missed:
###deterministically if their modeled height is less than whatever threshold;
###stochastically depending upon size. 
FalseGuys$y<-1
for (i in 1:nrow(TrueGuys)){
  TrueGuys$y[i]<-ifelse(TrueGuys$height[i]<.1, 0, rbinom(1, 1, plogis(-3+5*TrueGuys$height[i])))
}

###The above is a little tricky in that the coefficients are not easily amenable to a prior:
###(don't want some N(0, 10) prior for the effect of height on detectable, for example).
###Seems like this can be simplified.

###Just to give a sense of the association
#plot(TrueGuys$height, plogis(-3+15*TrueGuys$height)) 
#plot(TrueGuys$realHt, TrueGuys$y)
#plot(TrueGuys$height, TrueGuys$y)

###Could also add some simulated measure of confidence in the plant
###following what might be produced by a CNN.

###Some number of 'grid cells' are validated in some way in the field.
###Many ways to consider a sampling scheme for the validation.

###Anyway, here, let's imagine 200 cells randomly sampled.
Validation<-sample(grid$cell, 200)

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
Vals$z<-c(rep(0, sum(FalseGuys$cellID %in% Validation)), rep(1, sum(TrueGuys$cellID %in% Validation)))

###We don't see any of the NonVals if modeled height is <.1 m
NonVals<-NonVals[NonVals$height >=.1,]

##Can't falsely detect a Val if the modeled height is <.1 m either

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

###need to re-index the cells
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
           height2=Vals$height)
#validated=Validation, nonvalidated=grid$cell %notin% Validation) 

hinit<-c(rep(NA, nrow(NonVals)), rgamma(1000, 2, 10)) 
Zinit<-c(rep(1, nrow(NonVals)), rep(0, 1000))
Cinit<-Zinit+1
Inits <- list(z=Zinit, cl=Cinit, height=hinit, psi=0.5, p10=0.75, alpha_0=0, alpha_1=0, beta=matrix(0, 2, 2), shape=rep(5, 2), rate=rep(10, 2))

###need to check initial values...geting some bad logProbs
Shrub <- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                     data=Data, inits = Inits)

mcmcConf <- configureMCMC(Shrub, monitors = c("alpha_0", "alpha_1",  "N", "beta", 
                                              "psi", "p10", "shape", "rate"), useConjugacy=FALSE)

Rmcmc<-buildMCMC(mcmcConf)


compMCMC <- compileNimble(Rmcmc, Shrub)

samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=8000, nburnin=5000, thin=1, 
               nchains=1)

