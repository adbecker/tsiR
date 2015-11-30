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
#' @param fittype, the type of fit used. Options are 'all' which fits beta, sbar, and alpha, 
#' 'fixalpha', which fixes alpha at 0.97 and estimates beta and sbar, and
#' 'less' which fits only beta and fixes alpha at 0.97.
#' @param fit, the fitting method used. Options are 'glm' or 'bayesglm' which is a cheap adaptation
#' to include some bayesian approaches. For 'bayesglm' we use a gaussian prior with mean 1e-4.
#' @param sbar, the mean number of susceptibles. Only used if fittype='less'. Defaults to 0.05*mean(pop).
#' @param printon, whether to show diagnostic prints or not, defaults to FALSE.

estpars <- function(data, xreg = 'cumcases',IP = 2,
                    regtype = 'gaussian',sigmamax = 3,
                    userYhat = numeric(),fit='glm',sbar=0.05,
                    fittype = 'all',printon=F){
  
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
  
  Iadjusted <- data$cases * adj.rho
  
  datacopy <- data
  data$cases[data$cases ==0] <- 1
  
  period <- rep(1:(52/IP), round(nrow(data)+1))[1:(nrow(data)-1)]
  
  if(IP == 1){
    
    period <- rep(1:(52/2),each=2, round(nrow(data)+1))[1:(nrow(data)-1)]
    
  }
  
  Inew <- tail(Iadjusted,-1)+1
  lIold <- log(head(Iadjusted,-1)+1)
  Zold <- head(Z,-1)
  
  pop <- data$pop
  
  Smean <- seq(0.001, 0.4, by=0.001)*mean(pop)
  
  llik <- rep(NA, length(Smean))
  if(fittype == 'all'){
    
    for(i in 1:length(Smean)){
      lSold <- log(Smean[i] + Zold)
      
      if(fit == 'glm'){
        glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIold) + offset(lSold),
                      family=gaussian(link='log'))
      }
      
      if(fit == 'bayesglm'){
        glmfit <- bayesglm(Inew ~ -1 +as.factor(period) + (lIold) + offset(lSold),
                           family=gaussian(link='log'),prior.df=Inf,prior.mean=1e-4)
      }
      
      
      llik[i] <- glmfit$deviance
      
    }
    
    sbar <- Smean[which.min(llik)]
    
    lSold <- log(sbar + Zold)
    if(fit == 'glm'){
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ (lIold) + offset(lSold),
                    family=gaussian(link='log'))
    }
    
    if(fit == 'bayesglm'){
      
      glmfit <- bayesglm(Inew ~ -1 +as.factor(period)+ (lIold) + offset(lSold),
                         family=gaussian(link='log'),prior.df=Inf,prior.mean=1e-4)
    }
    
    beta <- exp(head(coef(glmfit),-1))
    alpha <- tail(coef(glmfit),1)
  }
  
  
  if(fittype == 'fixalpha'){
    
    alpha <- 0.97
    
    for(i in 1:length(Smean)){
      lSold <- log(Smean[i] + Zold)
      if(fit == 'glm'){
        
        
        glmfit <- glm(Inew ~ -1 +as.factor(period) + offset(alpha*lIold) + offset(lSold),
                      family=gaussian(link='log'))
        
      }
      
      if(fit == 'bayesglm'){
        
        
        glmfit <- bayesglm(Inew ~ -1 +as.factor(period) + offset(alpha*lIold) + offset(lSold),
                           family=gaussian(link='log'),prior.df=Inf,prior.mean=1e-4)
        
      }
      
      
      llik[i] <- glmfit$deviance
      
    }
    
    sbar <- Smean[which.min(llik)]
    
    lSold <- log(sbar + Zold)
    
    if(fit == 'glm'){
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIold) + offset(lSold),
                    family=gaussian(link='log'))
      
    }
    
    if(fit == 'bayesglm'){
      
      glmfit <- bayesglm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIold) + offset(lSold),
                         family=gaussian(link='log'),prior.df=Inf,prior.mean=1e-4)
      
    }
    
    
    beta <- exp(coef(glmfit))
  }
  
  
  if(fittype == 'less'){
    sbar <- sbar * mean(pop)
    alpha <- 0.97
    lSold <- log(sbar + Zold)
    
    if(fit == 'glm'){
      
      glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIold) + offset(lSold),
                    family=gaussian(link='log'))
      
    }
    
    if(fit == 'bayesglm'){
      
      glmfit <- bayesglm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIold) + offset(lSold),
                         family=gaussian(link='log'),prior.df=Inf,prior.mean=1e-4)
      
    }
    
    beta <- exp(coef(glmfit))
    
  }
  
  return(list('beta'=beta,'rho'=adj.rho,'Z'=Z,'sbar'=sbar,'alpha'=alpha,'loglik'=llik))
  
}