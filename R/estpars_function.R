#' estpars function
#' 
#' This function computes the set up to run the TSIR model, i.e. reconstructes susecptibles and 
#' estimates beta and alpha.
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
#' @param family, the family in the GLM regression, options are poisson and gaussian both with
#' log link. Default is Poisson.
#' @param fit, the fitting method used. Options are 'glm' or 'bayesglm' which is a cheap adaptation
#' to include some bayesian approaches. For 'bayesglm' we use a gaussian prior with mean 1e-4.
#' @param seasonality, the type of contact to use. Options are standard for 52/IP point contact or schoolterm for just a two point on off contact. Defaults to standard.
#' @param sbar, the mean number of susceptibles. Defaults to NULL, i.e. the function estimates sbar.
#' @param alpha, the mixing parameter. Defaults to NULL, i.e. the function estimates alpha.
#' @param printon, whether to show diagnostic prints or not, defaults to FALSE.

estpars <- function(data, xreg = 'cumcases',IP = 2,seasonality='standard',
                    regtype = 'gaussian',sigmamax = 3,family='gaussian',
                    userYhat = numeric(),alpha=NULL,sbar=NULL,
                    printon=F){
  
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
        Z <- residual.cases(Yhat,Y)
        rho <- derivative(X,Yhat)
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
  
  Inew <- tail(Iadjusted,-1)+1
  lIminus <- log(head(Iadjusted,-1)+1)
  Zminus <- head(Z,-1)
  
  pop <- data$pop
  
  minSmean <- max(0.01*pop,-(min(Z)+1))
  Smean <- seq(minSmean, 0.4*mean(pop), length=250)
  
  loglik <- rep(NA, length(Smean))
  
  if(length(input.alpha) == 0 && length(input.sbar) == 0){
    
    if(family == 'gaussian'){
      
      for(i in 1:length(Smean)){
        lSminus <- log(Smean[i] + Zminus)
        
        glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIminus) + offset(lSminus),
                      family=gaussian(link='log'))
        
        loglik[i] <- glmfit$deviance
        
      }
      
      sbar <- Smean[which.min(loglik)]
      
      lSminus <- log(sbar + Zminus)
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ (lIminus) + offset(lSminus),
                    family=gaussian(link='log'))
      
      
    }
    
    
    if(family == 'poisson'){
      
      Inew <- round(Inew)
      
      for(i in 1:length(Smean)){
        lSminus <- log(Smean[i] + Zminus)
        
        glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIminus) + offset(lSminus),
                      family=poisson(link='log'))
        
        loglik[i] <- glmfit$deviance
        
      }
      
      sbar <- Smean[which.min(loglik)]
      
      lSminus<- log(sbar + Zminus)
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ (lIminus) + offset(lSminus),
                    family=poisson(link='log'))
      
      
    }
    
    
    beta <- exp(head(coef(glmfit),-1))
    alpha <- tail(coef(glmfit),1)
  }
  
  
  if(length(input.alpha) == 1 && length(input.sbar) == 0){
    
    if(family == 'gaussian'){
      
      for(i in 1:length(Smean)){
        lSminus <- log(Smean[i] + Zminus)
        
        
        glmfit <- glm(Inew ~ -1 +as.factor(period) + offset(alpha*lIminus) + offset(lSminus),
                      family=gaussian(link='log'))
        
        
        
        
        loglik[i] <- glmfit$deviance
        
      }
      
      sbar <- Smean[which.min(loglik)]
      
      lSminus <- log(sbar + Zminus)
      
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIminus) + offset(lSminus),
                    family=gaussian(link='log'))
      
    }
    
    
    if(family == 'poisson'){
      
      Inew <- round(Inew)
      
      
      for(i in 1:length(Smean)){
        lSminus <- log(Smean[i] + Zminus)
        
        
        glmfit <- glm(Inew ~ -1 +as.factor(period) + offset(alpha*lIminus) + offset(lSminus),
                      family=poisson(link='log'))
        
        
        loglik[i] <- glmfit$deviance
        
      }
      
      sbar <- Smean[which.min(loglik)]
      
      lSminus <- log(sbar + Zminus)
      
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIminus) + offset(lSminus),
                    family=poisson(link='log'))
      
      
    }
    
    beta <- exp(coef(glmfit))
  }
  
  
  if(length(input.alpha) == 0 && length(input.sbar) == 1){
    
    sbar <- sbar * mean(pop)
    lSminus <- log(sbar + Zminus)
    
    if(family == 'gaussian'){
      
      
      glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIminus) + offset(lSminus),
                    family=gaussian(link='log'))
      
    }
    
    
    if(family == 'poisson'){
      
      Inew <- round(Inew)
      
      glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIminus) + offset(lSminus),
                    family=poisson(link='log'))
      
      
    }
    
    
    
    beta <- exp(head(coef(glmfit),-1))
    alpha <- tail(coef(glmfit),1)
  }
  
  
  if(length(input.alpha) == 1 && length(input.sbar) == 1){
    
    sbar <- sbar * mean(pop)
    lSminus <- log(sbar + Zminus)
    
    if(family == 'gaussian'){
      
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIminus) + offset(lSminus),
                    family=gaussian(link='log'))
      
      
    }
    
    if(family == 'poisson'){
      
      Inew <- round(Inew)
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIminus) + offset(lSminus),
                    family=poisson(link='log'))
      
      
    }
    
    beta <- exp(coef(glmfit))
    
  }
  
  
  return(list('X'=X,'Y'=Y,'Yhat'=Yhat,'Smean'=Smean,
              'beta'=head(beta[period],52/IP),'rho'=adj.rho,'Z'=Z,
              'sbar'=sbar,'alpha'=alpha,'loglik'=loglik))
  
}