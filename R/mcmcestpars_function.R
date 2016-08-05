#' @title mcmcestpars
#'
#' @description This function computes the set up to run the TSIR model, i.e. reconstructes susecptibles and
#' estimates beta and alpha using MCMC computations. Used the same way as estpars.
#'
#' @param data The data frame containing cases and interpolated births and populations.
#' @param xreg The x-axis for the regression. Options are 'cumcases' and 'cumbirths'. Defaults to 'cumcases'.
#' @param IP The infectious period in weeks. Defaults to 2 weeks.
#' @param regtype The type of regression used in susceptible reconstruction.
#' Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom),
#' 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector.
#' Defaults to 'gaussian' and if that fails then defaults to loess.
#' @param sigmamax The inverse kernal width for the gaussian regression. Default is 3.
#' Smaller, stochastic outbreaks tend to need a lower sigma.
#' @param userYhat The inputed regression vector if regtype='user'. Defaults to NULL.
#' @param n.chains Number of MCMC chains to use. Default is 3.
#' @param update.iter Number of MCMC iterations to use in the update aspect. Default is 10000.
#' @param n.iter Number of MCMC iterations to use. Default is 30000.
#' @param n.adapt Adaptive number for MCMC. Default is 1000.
#' @param burn.in Burn in number. Default is 100.
#' @param seasonality The type of contact to use. Options are standard for 52/IP point contact or schoolterm for just a two point on off contact or none for a single contact parameter. Defaults to standard.
#' @param sbar The mean number of susceptibles. Defaults to NULL, i.e. the function estimates sbar.
#' @param alpha The mixing parameter. Defaults to NULL, i.e. the function estimates alpha.
#' @param printon Whether to show diagnostic prints or not, defaults to FALSE.

mcmcestpars <- function(data, xreg = 'cumcases',IP = 2,
                        regtype = 'gaussian',sigmamax = 3,
                        seasonality='standard',
                        userYhat = numeric(),
                        update.iter=10000,
                        n.iter=30000, n.chains=3,
                        n.adapt=1000,burn.in=100,
                        sbar=NULL,alpha=NULL,
                        printon=F){


  rjagscheck <- 'rjags' %in% installed.packages()[,"Package"]
  if(rjagscheck == FALSE){
    stop('Package "rjags" is not installed, please install prior to using MCMC portions of code')
  }

  datacheck <- c('time','cases','pop','births')
  if(sum(datacheck %in% names(data)) < length(datacheck)){
    stop('data frame must contain "time", "cases", "pop", and "births" columns')
  }

  na.casescheck <- sum(is.na(data$cases))
 if(na.casescheck > 0){
   stop('there cannot be any NAs in the cases vector -- please correct')
 }

 na.birthscheck <- sum(is.na(data$births))
 if(na.casescheck > 0){
   stop('there cannot be any NAs in the births vector -- please correct')
 }

  xregcheck <- c('cumcases','cumbirths')
  if(xreg %in% xregcheck == F){
    stop('xreg must be either "cumcases" or "cumbirths"')
  }

  regtypecheck <- c('gaussian','lm','spline','lowess','loess','user')
  if(regtype %in% regtypecheck == F){
    stop("regtype must be one of 'gaussian','lm','spline','lowess','loess','user'")
  }

  if(length(sbar) == 1){
    if(sbar > 1 || sbar < 0){
      stop("sbar must be a percentage of the population, i.e. between zero and one.")
    }
  }

  seasonalitycheck <- c('standard','schoolterm','none')
  if(seasonality %in% seasonalitycheck == F){
    stop("seasonality must be either 'standard' or 'schoolterm' or 'none'")
  }

  if(n.iter < 5000){

    print('number of MCMC iterations less than 5000 -- increase')

  }

  print('MCMC may take a while')

  input.alpha <- alpha
  input.sbar <- sbar

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

  x <- seq(X[1], X[length(X)], length=length(X))
  y <- approxfun(X, Y)(x)
  y[1] <- y[2] - (y[3]-y[2])

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
        rho <- derivative(X,Yhat)
        Z <- residual.cases(Yhat,Y)
        if(length(which(rho<=1))==0){
          break()
        }
      }
      if(xreg == 'cumbirths'){
        rho <- derivative(X,Yhat)
        Z <- residual.births(rho,Yhat,Y)
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
    Z <- residual.births(rho,Yhat,Y)
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


  if(length(which(adj.rho < 1 )) > 1){
    warning('Reporting exceeds 100% -- use different regression')
  }

  Iadjusted <- data$cases * adj.rho

  datacopy <- data

  if(seasonality == 'standard'){

    period <- rep(1:(52/IP), round(nrow(data)+1))[1:(nrow(data)-1)]

    if(IP == 1){

      period <- rep(1:(52/2),each=2, round(nrow(data)+1))[1:(nrow(data)-1)]

    }

  }

  if(seasonality == 'schoolterm'){

    ## do school time in base two weeks and then interpolate
    term <- rep(1,26)
    term[c(1,8,15,16,17,18,19,23,26)] <- 2

    iterm <- round(approx(term,n=52/IP)$y)
    period <- rep(iterm, round(nrow(data)+1))[1:(nrow(data)-1)]

  }

  if(seasonality == 'none'){

    period <- rep(1,nrow(data)-1)
    period[nrow(data)-1] <- 2

  }

  Inew <- tail(Iadjusted,-1)+1
  lIminus <- log(head(Iadjusted,-1)+1)
  Zminus <- head(Z,-1)

  pop <- data$pop

  minSmean <- max(0.01*pop,-(min(Z)+1))
  Smean <- seq(minSmean, 0.4*mean(pop), length=250)

  alphalow <- NA
  alphahigh <- NA

  sbarlow <- NA
  sbarhigh <- NA

  ## even though dont profile to get sbar
  ## keep this in there for plotting function
  loglik <- rep(NA, length(Smean))
  if(length(input.alpha) == 0 && length(input.sbar) == 0){

    factorperiod <- as.factor(period)
    mod <- model.matrix(~-1+factorperiod)

    numseas <- length(unique(period))
    mymodel <- textConnection('model{
                              alpha ~ dunif(0.5,0.99)
                              for(season in 1:numseas){
                              beta[season] ~ dunif(-13,-3)
                              }
                              sbar ~ dunif(minSmean, 0.4*mean(pop))

                              sigma ~ dunif(0,10)

                              for (t in 1:N){

                              ## no intercept
                              regsum[t] <- mod[t,] %*%beta + alpha*lIminus[t] + log(Zminus[t]+sbar) + e[t]
                              rate[t] <- exp(regsum[t])
                              Inew[t] ~ dpois(rate[t])
                              e[t] ~ dnorm(0, (1/sigma^2))
                              }

  }')

    jags_data_list=list(
      "mod" = mod,
      "Inew"=round(Inew),
      "lIminus"=as.numeric(lIminus),
      "Zminus"=as.numeric(Zminus),
      "numseas" = numseas,
      'pop'=pop,
      'minSmean'=minSmean,
      "N" = length(lIminus)
    )

    theModel <- jags.model(mymodel,data=jags_data_list,n.chains=n.chains)
    update(theModel,update.iter)
    inits = list("alpha" = 0.97)
    mcmcsamples <- coda.samples(theModel,c("alpha","beta",'sigma','sbar'),
                                inits=inits,n.iter=n.iter, n.adapt=n.adapt,burn.in=burn.in)

    results <-  as.data.frame(mcmcsamples[[1]])

    mcmctruncated <- tail(mcmcsamples,5000)

    jagsres <- jagsfilter(mcmcresults = mcmctruncated)

    beta <- exp(jagsres[2:(length(unique(period))+1),1])
    betalow <- exp(jagsres[2:(length(unique(period))+1),2])
    betahigh <- exp(jagsres[2:(length(unique(period))+1),3])

    alpha <- jagsres[1,1]
    alphalow <- jagsres[1,2]
    alphahigh <- jagsres[1,3]

    sbar <- jagsres[length(unique(period))+2,1]
    sbarlow <- jagsres[length(unique(period))+2,2]
    sbarhigh <- jagsres[length(unique(period))+2,3]

  }

  if(length(input.alpha) == 1 && length(input.sbar) == 0){

    factorperiod <- as.factor(period)
    mod <- model.matrix(~-1+factorperiod)

    numseas <- length(unique(period))
    mymodel <- textConnection('model{
                            for(season in 1:numseas){
                            beta[season] ~ dunif(-13,-3)
                            }
                            sbar ~ dunif(minSmean, 0.4*mean(pop))

                            sigma ~ dunif(0,10)

                            for (t in 1:N){

                            ## no intercept
                            regsum[t] <- mod[t,] %*%beta + alpha*lIminus[t] + log(Zminus[t]+sbar) + e[t]
                            rate[t] <- exp(regsum[t])
                            Inew[t] ~ dpois(rate[t])
                            e[t] ~ dnorm(0, (1/sigma^2))
                            }

}')

    jags_data_list=list(
      "mod" = mod,
      "Inew"=round(Inew),
      "lIminus"=as.numeric(lIminus),
      "Zminus"=as.numeric(Zminus),
      "numseas" = numseas,
      'alpha' = alpha,
      'pop'=pop,
      'minSmean'=minSmean,
      "N" = length(lIminus)
    )

    theModel <- jags.model(mymodel,data=jags_data_list,n.chains=n.chains)
    update(theModel,update.iter)
    mcmcsamples <- coda.samples(theModel,c("beta",'sigma','sbar'),
                                inits=inits,n.iter=n.iter, n.adapt=n.adapt,burn.in=burn.in)

    results <-  as.data.frame(mcmcsamples[[1]])

    mcmctruncated <- tail(mcmcsamples,5000)

    jagsres <- jagsfilter(mcmcresults = mcmctruncated)

    beta <- exp(jagsres[1:(length(unique(period))),1])
    betalow <- exp(jagsres[1:(length(unique(period))),2])
    betahigh <- exp(jagsres[1:(length(unique(period))),3])

    sbar <- jagsres[length(unique(period))+1,1]
    sbarlow <- jagsres[length(unique(period))+1,2]
    sbarhigh <- jagsres[length(unique(period))+1,3]

  }

  if(length(input.alpha) == 0 && length(input.sbar) == 1){

    sbar <- sbar*mean(pop)

    factorperiod <- as.factor(period)
    mod <- model.matrix(~-1+factorperiod)

    numseas <- length(unique(period))
    mymodel <- textConnection('model{
                            for(season in 1:numseas){
                            beta[season] ~ dunif(-13,-3)
                            }
                            alpha ~ dunif(0.3, 0.99)

                            sigma ~ dunif(0,10)

                            for (t in 1:N){

                            ## no intercept
                            regsum[t] <- mod[t,] %*%beta + alpha*lIminus[t] + log(Zminus[t]+sbar) + e[t]
                            rate[t] <- exp(regsum[t])
                            Inew[t] ~ dpois(rate[t])
                            e[t] ~ dnorm(0, (1/sigma^2))
                            }

}')

    jags_data_list=list(
      "mod" = mod,
      "Inew"=round(Inew),
      "lIminus"=as.numeric(lIminus),
      "Zminus"=as.numeric(Zminus),
      "numseas" = numseas,
      'sbar' = sbar,
      "N" = length(lIminus)
    )

    theModel <- jags.model(mymodel,data=jags_data_list,n.chains=n.chains)
    update(theModel,update.iter)
    mcmcsamples <- coda.samples(theModel,c("beta",'sigma','alpha'),
                                inits=inits,n.iter=n.iter, n.adapt=n.adapt,burn.in=burn.in)

    results <-  as.data.frame(mcmcsamples[[1]])

    mcmctruncated <- tail(mcmcsamples,5000)

    jagsres <- jagsfilter(mcmcresults = mcmctruncated)

    beta <- exp(jagsres[2:(length(unique(period))+1),1])
    betalow <- exp(jagsres[2:(length(unique(period))+1),2])
    betahigh <- exp(jagsres[2:(length(unique(period))+1),3])

    alpha <- jagsres[1,1]
    alphalow <- jagsres[1,2]
    alphahigh <- jagsres[1,3]

  }

  if(length(input.alpha) == 1 && length(input.sbar) == 1){

    sbar <- sbar * mean(pop)
    lSminus <- log(sbar + Zminus)

    lSminus[is.nan(lSminus)] <- 0
    lSminus[lSminus < 0] <- 0

    factorperiod <- as.factor(period)
    mod <- model.matrix(~-1+factorperiod)

    numseas <- length(unique(period))
    mymodel <- textConnection('model{
                            for(season in 1:numseas){
                            beta[season] ~ dunif(-13,-3)
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

    jagsres <- jagsfilter(mcmcresults = mcmctruncated)

    beta <- exp(jagsres[1:length(unique(period)),1])
    betalow <- exp(jagsres[1:length(unique(period)),2])
    betahigh <- exp(jagsres[1:length(unique(period)),3])

  }

  if(seasonality == 'none'){
    beta[2] <- beta[1]
    beta <- mean(beta)
    period <- rep(1,nrow(data)-1)
  }

  contact <- as.data.frame(cbind('time'=seq(1,length(beta[period]),1),
                                 betalow[period],beta[period],betahigh[period]),row.names=F)
  names(contact) <- c('time','betalow','beta','betahigh')
  contact <- head(contact,52/IP)
 

  return(list('X'=X,'Y'=Y,'Yhat'=Yhat,
              'mcmcsamples'=mcmcsamples,
              'beta'=contact$beta,'contact'=contact,'rho'=adj.rho,'pop'=pop,
              'Z'=Z,'sbar'=sbar,'sbarlow'=sbarlow,'sbarhigh'=sbarhigh,
              'alpha'=alpha, 'alphalow'=alphalow,'alphahigh'=alphahigh,
              'loglik'=loglik,'Smean'=Smean
  ))

}
