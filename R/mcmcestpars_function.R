#' mcmcestpars function
#' 
#' This function computes the set up to run the TSIR model, i.e. reconstructes susecptibles and 
#' estimates beta and alpha using MCMC computations.
#' 
#' @param data, the data frame containing cases and interpolated births and populations.
#' @param xreg, the x-axis for the regression. Options are 'cumcases' and 'cumbirths'. Defaults to 'cumcases'.
#' @param IP, the infectious period in weeks. Defaults to 2 weeks.
#' @param regtype, the type of regression used in susceptible reconstruction. 
#' Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom),
#' 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector.
#' Defaults to 'gaussian' and if that fails then defaults to loess.
#' @param sigmamax, the inverse kernal width for the gaussian regression. Default is 3. 
#' Smaller, stochastic outbreaks tend to need a lower sigma.
#' @param userYhat, the inputed regression vector if regtype='user'. Defaults to NULL.
#' @param fittype, the type of fit used. Options are 'all' which fits beta, sbar, and alpha, 
#' 'fixalpha', which fixes alpha at 0.97 and estimates beta and sbar, and
#' 'less' which fits only beta and fixes alpha at 0.97.
#' @param n.chains, number of MCMC chains to use. Default is 3.
#' @param update.iter, number of MCMC iterations to use in the update aspect. Default is 10000.
#' @param n.iter, number of MCMC iterations to use. Default is 30000.
#' @param n.adapt, adaptive number for MCMC. Default is 1000.
#' @param burn.in, burn in number. Default is 100.
#' @param sbar, the mean number of susceptibles. Only used if fittype='less'. Defaults to 0.05*mean(pop).
#' @param printon, whether to show diagnostic prints or not, defaults to FALSE.

mcmcestpars <- function(data, xreg = 'cumcases',IP = 2,
                        regtype = 'gaussian',sigmamax = 3,
                        userYhat = numeric(),
                        update.iter=10000,
                        n.iter=30000, n.chains=3, 
                        n.adapt=1000,burn.in=100,
                        sbar=0.05,
                        fittype = 'all',printon=F){
  
  if(n.iter < 5000){
    
    print('number of MCMC iterations less than 5000 -- increase')
    
  }
  
  print('MCMC may take a while')
  
  cumbirths <- cumsum(data$births)
  cumcases <- cumsum(data$cases)
  
  if(xreg == 'cumcases'){
    X <- cumcases
    Y <- cumbirths
  }
  
  if(xreg == 'cumbirths'){
    X <- cumbirths
    Y <- cumcases
  }
  
  x = linspace(X[1], X[length(X)], length(X))
  y = approxfun(X, Y)(x)
  y[1] = y[2] - (y[3]-y[2])
  
  if(regtype == 'lm'){
    Yhat <- predict(lm(Y~X))
  }
  
  if(regtype == 'lowess'){
    Yhat <- lowess(X,Y,f = 2/3, iter = 1)$y
  }
  
  if(regtype == 'loess'){
    Yhat <- predict(loess(y~x,se=T,family='gaussian',degree=1,model=T),X)
  }
  
  if(regtype == 'spline'){
    Yhat <- predict(smooth.spline(x,y,df=2.5),X)$y
  }
  
  if(regtype == 'gaussian'){
    
    sigvec <- seq(sigmamax,0,-0.1)
    for(it in 1:length(sigvec)){
      
      if(printon == T){
        print(sprintf('gaussian regression attempt number %d',it))
      }
      
      Yhat <- predict(gausspr(x,y,variance.model=T,fit=T,tol=1e-7,
                              var=9.999999999999999999e-3,
                              kernel="rbfdot",
                              kpar=list(sigma=sigvec[it])),X)
      
      
      if(sigvec[it] <= min(sigvec)){
        ## use the loess then
        print('guassian regressian failed -- switching to loess regression')
        Yhat <- predict(loess(y~x,se=T,family='gaussian',degree=1,model=T),X)
      }
      
      
      if(xreg == 'cumcases'){
        Z <- residual.cases(Yhat,Y)
        rho <- derivative(X,Yhat)
        if(length(which(rho<=1))==0){
          break()
        }
      }
      if(xreg == 'cumbirths'){
        Z <- residual.births(Yhat,Y)
        rho <- derivative(X,Yhat)
        if(length(which(rho>=1))==0 && length(which(rho<0)) == 0){
          break()
        }
      }
    }
  }
  
  
  if(regtype == 'user'){
    Yhat <- userYhat
    if(length(Yhat)==0){
      stop('Yhat returns numeric(0) -- make sure to input a userYhat under regtype=user')
    }
  }
  
  rho <- derivative(X,Yhat)
  
  if(xreg == 'cumcases'){
    Z <- residual.cases(Yhat,Y)
  }
  
  if(xreg == 'cumbirths'){
    Z <- residual.births(Yhat,Y)
  }
  
  if(xreg == 'cumcases'){
    adj.rho <- rho
  }
  if(xreg == 'cumbirths'){
    adj.rho <- 1/rho
  }
  
  if(regtype == 'lm'){
    adj.rho <- signif(adj.rho,3)
    
  }
  
  if(adj.rho < 1){
    stop()
  }
  
  if(length(which(adj.rho < 1 )) > 1){
    stop('Reporting exceeds 100% -- use different regression')
  }
  
  Iadjusted <- data$cases * adj.rho
  
  datacopy <- data
  
  period <- rep(1:(52/IP), round(nrow(data)+1))[1:(nrow(data)-1)]
  
  if(IP == 1){
    
    period <- rep(1:(52/2),each=2, round(nrow(data)+1))[1:(nrow(data)-1)]
    
  }
  
  Inew <- tail(Iadjusted,-1)+1
  lIminus <- log(head(Iadjusted,-1)+1)
  Zminus <- head(Z,-1)
  
  pop <- data$pop
  
  Smean <- seq(0.001, 0.4, by=0.001)*mean(pop)
  
  alphalow <- NA
  alphahigh <- NA
  
  loglik <- rep(NA, length(Smean))
  if(fittype == 'all'){
    
    
    for(i in 1:length(Smean)){
      lSminus <- log(Smean[i] + Zminus)
      
      glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIminus) + offset(lSminus),
                    family=poisson(link='log'))
      
      
      loglik[i] <- glmfit$deviance
      
    }
    
    sbar <- Smean[which.min(loglik)]
    
    lSminus <- log(sbar + Zminus)
    
    lSminus[is.nan(lSminus)] <- 0
    lSminus[lSminus < 0] <- 0
    
    
    factorperiod <- as.factor(period)
    mod <- model.matrix(~-1+factorperiod)
    
    numseas <- length(unique(period))
    mymodel <- textConnection('model{
                                alpha ~ dunif(0.5,0.99)
                                for(season in 1:numseas){
                                beta[season] ~ dunif(-12,-3)
                                }
                                
                                sigma ~ dunif(0,10)
                                
                                for (t in 1:N){
                                
                                ## no intercept
                                regsum[t] <- mod[t,] %*%beta + alpha*lIminus[t] + lSminus[t] + e[t]
                                rate[t] <- exp(regsum[t])
                                Inew[t] ~ dpois(rate[t])
                                e[t] ~ dnorm(0, (1/sigma^2))
                                }
                                
    }')
    
    jags_data_list=list(
      "mod" = mod,
      "Inew"=round(Inew),
      "lIminus"=as.numeric(lIminus),
      "lSminus"=as.numeric(lSminus),
      "numseas" = numseas,
      "N" = length(lIminus)
    )
    
    theModel <- jags.model(mymodel,data=jags_data_list,n.chains=n.chains)
    update(theModel,update.iter)
    inits = list("alpha" = 0.97)
    mcmcsamples <- coda.samples(theModel,c("alpha","beta",'sigma'),
                                inits=inits,n.iter=n.iter, n.adapt=n.adapt,burn.in=burn.in)
    
    results <-  as.data.frame(mcmcsamples[[1]])
    
    mcmctruncated <- tail(mcmcsamples,5000)
    
    jagsres <- jagsresults(x=mcmctruncated, param=names(results))
    jagsres <- jagsres[,c('mean', '2.5%', '97.5%')]
    
    beta <- exp(jagsres[2:(length(unique(period))+1),1])
    betalow <- exp(jagsres[2:(length(unique(period))+1),2])
    betahigh <- exp(jagsres[2:(length(unique(period))+1),3])
    
    alpha <- jagsres[1,1]
    alphalow <- jagsres[1,2]
    alphahigh <- jagsres[1,3]
  }
  
  
  
  
  if(fittype == 'fixalpha'){
    
    alpha <- 0.97
    
    
    for(i in 1:length(Smean)){
      lSminus <- log(Smean[i] + Zminus)
      
      
      glmfit <- glm(Inew ~ -1 +as.factor(period) + offset(alpha*lIminus) + offset(lSminus),
                    family=poisson(link='log'))
      
      
      loglik[i] <- glmfit$deviance
      
    }
    
    
    
    sbar <- Smean[which.min(loglik)]
    
    lSminus <- log(sbar + Zminus)
    
    
    lSminus[is.nan(lSminus)] <- 0
    lSminus[lSminus < 0] <- 0
    
    
    factorperiod <- as.factor(period)
    mod <- model.matrix(~-1+factorperiod)
    
    numseas <- length(unique(period))
    mymodel <- textConnection('model{
                              for(season in 1:numseas){
                              beta[season] ~ dunif(-12,-3)
                              }
                              
                              sigma ~ dunif(0,10)
                              
                              for (t in 1:N){
                              
                              ## no intercept
                              regsum[t] <- mod[t,] %*%beta + alpha*lIminus[t] + lSminus[t] + e[t]
                              rate[t] <- exp(regsum[t])
                              Inew[t] ~ dpois(rate[t])
                              e[t] ~ dnorm(0, (1/sigma^2))
                              }
                              
  }')
    
    jags_data_list=list(
      "mod" = mod,
      "alpha" = alpha,
      "Inew"=round(Inew),
      "lIminus"=as.numeric(lIminus),
      "lSminus"=as.numeric(lSminus),
      "numseas" = numseas,
      "N" = length(lIminus)
    )
    
    theModel <- jags.model(mymodel,data=jags_data_list,n.chains=n.chains)
    update(theModel,update.iter)
    mcmcsamples <- coda.samples(theModel,c("beta",'sigma'),
                                n.iter=n.iter, n.adapt=n.adapt,burn.in=burn.in)
    
    results <-  as.data.frame(mcmcsamples[[1]])
    
    mcmctruncated <- tail(mcmcsamples,5000)
    
    jagsres <- jagsresults(x=mcmctruncated, param=names(results))
    jagsres <- jagsres[,c('mean', '2.5%', '97.5%')]
    
    beta <- exp(jagsres[1:length(unique(period)),1])
    betalow <- exp(jagsres[1:length(unique(period)),2])
    betahigh <- exp(jagsres[1:length(unique(period)),3])
    
    
  }
  
  
  
  if(fittype == 'less'){
    sbar <- sbar * mean(pop)
    alpha <- 0.97
    lSminus <- log(sbar + Zminus)
    
    
    
    lSminus[is.nan(lSminus)] <- 0
    lSminus[lSminus < 0] <- 0
    
    
    factorperiod <- as.factor(period)
    mod <- model.matrix(~-1+factorperiod)
    
    numseas <- length(unique(period))
    mymodel <- textConnection('model{
                              for(season in 1:numseas){
                              beta[season] ~ dunif(-12,-3)
                              }
                              
                              sigma ~ dunif(0,10)
                              
                              for (t in 1:N){
                              
                              ## no intercept
                              regsum[t] <- mod[t,] %*%beta + alpha*lIminus[t] + lSminus[t] + e[t]
                              rate[t] <- exp(regsum[t])
                              Inew[t] ~ dpois(rate[t])
                              e[t] ~ dnorm(0, (1/sigma^2))
                              }
                              
  }')
    
    jags_data_list=list(
      "mod" = mod,
      "alpha" = alpha,
      "Inew"=round(Inew),
      "lIminus"=as.numeric(lIminus),
      "lSminus"=as.numeric(lSminus),
      "numseas" = numseas,
      "N" = length(lIminus)
    )
    
    theModel <- jags.model(mymodel,data=jags_data_list,n.chains=n.chains)
    update(theModel,update.iter)
    mcmcsamples <- coda.samples(theModel,c("beta",'sigma'),
                                n.iter=n.iter, n.adapt=n.adapt,burn.in=burn.in)
    
    results <-  as.data.frame(mcmcsamples[[1]])
    
    mcmctruncated <- tail(mcmcsamples,5000)
    
    jagsres <- jagsresults(x=mcmctruncated, param=names(results))
    jagsres <- jagsres[,c('mean', '2.5%', '97.5%')]
    
    beta <- exp(jagsres[1:length(unique(period)),1])
    betalow <- exp(jagsres[1:length(unique(period)),2])
    betahigh <- exp(jagsres[1:length(unique(period)),3])
    
    
  }
  
  contact <- as.data.frame(cbind('time'=seq(1,length(beta),1),betalow,beta,betahigh))
  
  
  
  return(list('mcmcsamples'=mcmcsamples,
              'beta'=beta,'contact'=contact,'rho'=adj.rho,'pop'=pop,
              'Z'=Z,'sbar'=sbar,'alpha'=alpha,
              'alphalow'=alphalow,'alphahigh'=alphahigh,
              'res'=res,'loglik'=loglik
              ))
  
}