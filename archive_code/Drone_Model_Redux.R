###add covariance between size and density..
# SimDrone<-function(nsims, ncells, alphaTrue, betaTrue, alphaFalse, betaFalse1,
#                    betaFalse2, deltaFalse, gammaFalse,
#                    deltaTrue, gammaTrue, p, nVal){
#   N_True_Val<-rep(NA, nsims)
#   N_False_Val<-rep(NA, nsims)
#   N_True<-rep(NA, nsims)
#   N_est<-rep(NA, nsims)
#   beta_est<-matrix(NA, nsims, 2)
#   for (s in 1:nsims){
    #print(s)
    ###create landscape
    ncells=1600
    alphaTrue=-.5 ###intercept--expecting about .5 real plants per "average" cell--can modify!
    betaTrue=.7
    alphaFalse=0 ###intercept--expecting about 1 false plant per average cell--can modify!
    betaFalse1=1
    betaFalse2=-.5 
    deltaFalse=-.1 ##intercept for height. note, am assuming later that height is N(eta, sigma)
    gammaFalse = .5
    deltaTrue=0 ##intercept for height. 
    gammaTrue=1
    sigmaFalse=0.25
    sigmaTrue=0.25
    p=.9
    nVal=50 ###cells
    
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
    N_True<-sum(grid$TrueN)
    
    grid$FalseDens<-exp(alphaFalse+betaFalse1*grid$Cov+betaFalse2*grid$Cov^2)
    grid$FalseN<-rpois(nrow(grid), grid$FalseDens)
    
    ###"distribute" the individuals--cast the N into an individual data frame
    FalseGuys<-data.frame(ID=1:sum(grid$FalseN), cellID=rep(grid$cell, grid$FalseN))
    TrueGuys<-data.frame(ID=1:sum(grid$TrueN), cellID=rep(grid$cell, grid$TrueN))
    
    ###generate heights
    ###tricky in that the height will vary by state and cell...
    FalseGuys$height<-rnorm(nrow(FalseGuys), mean = deltaFalse+gammaFalse*grid$Cov[FalseGuys$cellID], sd=.25)
    TrueGuys$height<-rnorm(nrow(TrueGuys), mean = deltaTrue+gammaTrue*grid$Cov[TrueGuys$cellID], sd=.25)
    
    ###simulate detections
    FalseGuys$y<-1
    for (i in 1:nrow(TrueGuys)){
      TrueGuys$y[i]<-rbinom(1, 1, p)
    }
    
    ###Validation process...
    Validation<-sample(grid$cell, nVal)
    
    N_True_Val<-sum(TrueGuys$cellID %in% Validation)
    ###number of confirmed false guys
    N_False_Val<-sum(FalseGuys$cellID %in% Validation)
    ###Number of burned cells validated
    
    `%notin%` <- Negate(`%in%`)
    NonVals<-rbind(FalseGuys[FalseGuys$cellID %notin% Validation, ], TrueGuys[TrueGuys$cellID %notin% Validation, ])
    Vals<-rbind(FalseGuys[FalseGuys$cellID %in% Validation, ], TrueGuys[TrueGuys$cellID %in% Validation, ])
    
    NonVals$z<-NA
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
    
    ###all of the above is now, I think somewhat extraneous.
    ###easier to enter this as a number by cell with state[cell, n[cell]] and height[cell, n[cell]] where 
    ###n[cell]>0. Nimble sucks with logical statements, so have to do a bunch of subsetting/slicing here.
    
    ###number of things seen in non validated cells
    n<-rep(0, nrow(grid)-length(Validation))
    for (i in 1:length(n)){
      n[i]<-length(which(NonVals$s==i))
    }
    
    zerosNV<-which(n==0)
    cNV<-which(n>0)
    z<-matrix(0, length(cNV), max(n)) ##are the non-validated entities real or not
    height<-matrix(0, length(cNV), max(n)) ##what's the height of the non-validated entities
    for (i in 1:nrow(height)){
        height[i, 1:n[cNV[i]]]<-NonVals$height[which(NonVals$s==cNV[i])]
        z[i, 1:n[cNV[i]]]<-NA
      }
   
    nFP<-rep(0, length(Validation)) ###number of false things per validated cell
    nTP<-rep(0, length(Validation))  ###number of true things
    nV<-rep(0, length(Validation)) ###number of things
    for (i in 1:length(nFP)){
      nFP[i]<-length(which(Vals$s==i & Vals$z==0))
      nTP[i]<-length(which(Vals$s==i & Vals$z==1))
      nV[i]<-length(which(Vals$s==i))
    }
    zeroesV<-which(nFP+nTP==0)   ###which validated cells have nothing.
    cV<-which(nFP+nTP>0)   ##which validated cells have something
    z2<-matrix(NA, length(cV), max(nFP+nTP)) ###state, height, of the validated things
    height2<-matrix(0, length(cV), max(nFP+nTP))
    y<-matrix(0, length(cV), max(nFP+nTP)) ###were the validated things detected
    for (i in 1:nrow(height2)){
      height2[i, 1:nV[cV[i]]]<-Vals$height[which(Vals$s==cV[i])]
      z2[i, 1:nV[cV[i]]]<-Vals$z[which(Vals$s==cV[i])]
      y[i, 1:nV[cV[i]]]<-Vals$y[which(Vals$s==cV[i])] 
    }
    
    
    library(nimble)
    ###EM things for deriving latent N. Maybe not working as desired....
    GetN<-nimbleFunction(
      run=function(lam=double(0), p=double(0), n=integer(0)){
      returnType(integer(0))
      K<-round(5*sqrt(lam)+lam-n) 
      n0 <- 0:K ###potential number of missed individuals
      tmp1<-dpois(n+n0, lam, log=1) ###pr N= nseen+range of nmissed|E(n)
      tmp2<-dbinom(n, n+n0,p, log=1) ###pr seeing n given N=n+nmissed and p
      tmp <- exp(tmp1+tmp2)
      post<- tmp/sum(tmp)
      nprime<-rcat(n = 1, post) 
      N<-nprime+n
      return(N)
      })
    
    Mod<-nimbleCode({
      
      for (c in 1:2){
        for (b in 1:nfx[c]){
          beta[b, c]~dnorm(0, sd=1) ###describe intensity of point processes
        }
      }
      
      p~dunif(0, 1) ###pr real shrub detected from air
      for (c in 1:2){ ###relation between height and cov
        delta[c]~dnorm(0, sd=1)
        gamma[c]~dnorm(0, sd=1)
        sigma[c]~T(dnorm(0, .5), 0, )
      }
      
      ###nimble does not accept vectors of indices...quite tedious here.
      
      for (j in 1:nNaive){
        mu[j, 2]<-exp(inprod(beta[1:nfx[2], 2], XNaive[j, 1:nfx[2]]))
        ###difference between naive and val--only detected real guys in naive
        #probs[j, 2]<-mu[j, 2]/sum(mu[1:nNaive, 2])
        mu[j, 1]<-exp(inprod(beta[1:nfx[1], 1], XNaive[j, 1:nfx[1]]))  
        #probs[j, 1]<-mu[j, 1]/sum(mu[1:nNaive, 2])
        n[j]~dpois(p*mu[j, 2]+mu[j, 1]) ###nice thing about poisson--sums and thinning still marginally poisson
        psi[j]<-(p*mu[j, 2])/(mu[j, 2]+mu[j, 1])
        ###psi= probability individual is real in a given cell. 
        ### this may not be correct--maybe missing size....
      }
      
      for (j in 1:nVal){ ##counts in all validated cells
        mu2[j, 2]<-exp(inprod(beta[1:nfx[2], 2], XVal[j, 1:nfx[2]]))
        nTP[j]~dpois(mu2[j, 2])
        #probs2[j, 2]<-mu2[j, 2]/sum(mu2[1:nVal, 2])
        mu2[j, 1]<-exp(inprod(beta[1:nfx[1], 1], XVal[j, 1:nfx[1]]))
        nFP[j]~dpois(mu2[j, 1])
        #probs2[j, 1]<-mu2[j, 1]/sum(mu2[1:nVal, 1])
      }
      
      ###validated plants...
      for (j in 1:nonZeroVal){
        for (c in nV[cV[j]]){
          height2[j, c]~dnorm(delta[class2[j, c]]+gamma[class2[j, c]]*XVal[cV[j], 2], sd=sigma[class2[j, c]])
          y[j, c]~dbern(z2[j, c]*p+(1-z2[j, c]))
        }
      }
      
      for (j in 1:nonZeroNonVal){
        for (c in n[cNV[j]]){
          z[j, c]~dbern(psi[cNV[j]])
          ###related to comments on psi above--this is...maybe not right and not 
          ###capturing the size part correctly. Certainly, tracing nseen below
          ###suggests extremely little variability, which results in poor coverage
          cl[j, c]<-z[j, c]+1
          height[j, c]~dnorm(delta[cl[j, c]]+gamma[cl[j, c]]*XNaive[cNV[j], 2], sd=sigma[cl[j,c]])
        }
        nseen[j]<-sum(z[j, 1:maxn]) ###guess at the number of real individuals in the cell.
        ###We could derive this and z after the fact using EM (pr z[]==1|size, cell, etc.)
        ###Maybe easier.
    }
      ##All this is earlier junk
      # EN[1, 1]<-sum(mu[1:nNaive, 1]) ###Expected  abundance of fakes
      # EN[2, 1]<-p11*sum(mu[1:nNaive, 2]) ###Expected abundance of reals in unvalidated cells
      # EN[1, 2]<-sum(mu2[1:nVal, 1]) ###Expected abundance of fakes in validated cells
      # EN[2, 2]<-sum(mu2[1:nVal, 2]) ###expected abundance of reals in validated cells
      # psi[1]<-EN[2, 1]/(EN[2, 1]+EN[1, 1]) ##
      # psi[2]<-EN[2, 2]/(EN[2, 2]+EN[1, 2]) ##
      # 
      # 
      # 
      # #EN[2, 1]<-sum(mu[1:nNaive, 2]) ###Expected abundance of reals in unvalidated cells
      # #EN[2, 2]<-sum(mu2[1:nVal, 2]) ###expected abundance of reals in validated cells
      # #psi[1]<-sum(EN[2, 1])/M
      # #psi[2]<-sum(EN[2, 2])/V
      # #for (i in 1:M){ ###M this is the unverified data
      # #  z[i]~dbern(psi[1]) ###does the shrub exist
      # #  cl[i]<-z[i]+1
      # #  s[i]~dcat(probs[1:nNaive, cl[i]]) ###where--which pixel--does the shrub exist
      # #  height[i]~dnorm(delta[cl[i]]+gamma[cl[i]]*XNaive[s[i], 2], sd=sigma[cl[i]])
      # #  #y[i]~dbern(z[i]*p11+(1-z[i])*p10)
      # #}
      # for (i in 1:V){ ###M this is the verified data
      #   z2[i]~dbern(psi[2]) ###does the shrub exist
      #   class2[i]<-z2[i]+1
      #   s2[i]~dcat(probs2[1:nVal, class2[i]]) ###where--which pixel--does the shrub exist
      #   height2[i]~dnorm(delta[class2[i]]+gamma[class2[i]]*XVal[s2[i], 2], sd=sigma[cl[i]])
      #   y2[i]~dbern(z2[i]*p11+(1-z2[i]))
      # }
      N_seen<-sum(nseen[1:nonZeroNonVal]) ###guess at the number of unvalided shrubs that exist
      lam<-sum(mu[1:nNaive, 2])
      ###how many of the observed plants exist...
      ###if N[cell] is Poisson, sum(N[cell is poisson]), and so could do 
      ###some sort of joint poisson for P(n+nmissed|lambda) and NB(nmissed|n, p11)
      
    })
    
    
    Constants<-list(M=nrow(NonVals), V=nrow(Vals), 
                    nfx=c(3, 2), nNaive=nrow(grid)-length(Validation), nVal=length(Validation),
                    class2=z2+1, nonZeroVal=length(cV), nonZeroNonVal=length(cNV), maxn=max(n), nV=nV,
                    cV=cV, cNV=cNV)
    
    Data<-list(y=y, 
               z2=z2, height=height, height2=height2, 
               n=n, nFP=nFP, nTP=nTP, 
               XNaive=cbind(rep(1, nrow(grid)-length(Validation)), grid$Cov[-Validation], grid$Cov[-Validation]^2),
               XVal=cbind(rep(1, length(Validation)), grid$Cov[Validation], grid$Cov[Validation]^2)) 
    Zinit<-matrix(rbinom(dim(z)[1]*dim(z)[2], 1, .5), dim(z)[1], dim(z)[2])
    Zinit[z==0]<-0
    Cinit<-Zinit+1
    Inits <- list(z=Zinit, cl=Cinit, p=.75, beta=matrix(0, 3, 2), gamma=rep(0, 2),
                  delta=rep(0, 2))
    
    Shrub <- nimbleModel(code = Mod, name = 'Shrub', constants = Constants,
                         data=Data, inits = Inits)
    
    mcmcConf <- configureMCMC(Shrub, monitors = c("p", "beta", 
                                                   "delta", "gamma", "sigma", 
                                                   "N_seen", #"nseen", 
                                                   "lam"), useConjugacy=FALSE)
    
    Rmcmc<-buildMCMC(mcmcConf)
    
    
    compMCMC <- compileNimble(Rmcmc, Shrub)
    
    samps<-runMCMC(mcmc = compMCMC$Rmcmc,
                   niter=10000, nburnin=5000, thin=5, 
                   nchains=1)
    
    Nhats<-rep(NA, nrow(samps))
    
    for (s in 1:length(Nhats)){
      Nhats[s]<-GetN(samps[s, "lam"], samps[s, "p"], samps[s, "N_seen"])
    }
    hist(Nhats+sum(nTP))
    abline(v=sum(grid$TrueN), lwd=2)    ###this is,in theory, the 

    
    ##extract info from model...
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    N_est<-Mode(Nhats)+sum(nTP)
    ###probably should use 90% CI unless running for a long time....
    N_lci[s]<-quantile(Nhats, probs=.05)+sum(nTP)
    N_uci[s]<-quantile(Nhats, probs=.95)+sum(nTP)

