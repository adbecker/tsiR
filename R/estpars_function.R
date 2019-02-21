#' @title estpars
#' @description This function computes the set up to run the TSIR model, i.e. reconstructs susecptibles and
#' estimates beta and alpha. This can be plugged into simulatetsir.
#' @param data The data frame containing cases and interpolated births and populations.
#' @param xreg The x-axis for the regression. Options are 'cumcases' and 'cumbirths'. Defaults to 'cumcases'.
#' @param IP The infectious period in weeks. This should be the same as your timestep. Defaults to 2 weeks.
#' @param regtype The type of regression used in susceptible reconstruction.
#' Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom),
#' 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector.
#' Defaults to 'gaussian' and if that fails then defaults to loess.
#' @param sigmamax The inverse kernal width for the gaussian regression. Default is 3.
#' Smaller, stochastic outbreaks tend to need a lower sigma.
#' @param userYhat The inputed regression vector if regtype='user'. Defaults to NULL.
#' @param family The family in the GLM regression. One can use any of the GLM ones, but the options are essentially
#' 'poisson' (with link='log'), 'gaussian' (with link='log' or 'identity'), or 'quasipoisson' (with link='log'). Default is 'gaussian'.
#' @param link The link function used with the glm family. Options are link='log' or 'identity'. Default is 'identity'.
#' to include some bayesian approaches. For 'bayesglm' we use a gaussian prior with mean 1e-4.
#' @param seasonality The type of contact to use. Options are standard for 52/IP point contact or schoolterm for just a two point on off contact or none for a single contact parameter. Defaults to standard.
#' @param sbar The mean number of susceptibles. Defaults to NULL, i.e. the function estimates sbar.
#' @param alpha The mixing parameter. Defaults to NULL, i.e. the function estimates alpha.
#' @param printon Whether to show diagnostic prints or not, defaults to FALSE.
#' @examples
#' \dontrun{
#' require(kernlab)
#' London <- twentymeas[["London"]]
#' parms <- estpars(London)
#' names(parms)
#' sim <- simulatetsir(London,parms=parms,inits.fit=FALSE)
#'plotres(sim)
#'}

estpars <- function(data, xreg = 'cumcases',IP = 2,seasonality='standard',
                    regtype = 'gaussian',sigmamax = 3,family='gaussian',link='identity',
                    userYhat = numeric(),alpha=NULL,sbar=NULL,
                    printon=F){

  ## for better annotations please see runtsir function

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

  linkcheck <- c('log','identity')
  if(link %in% linkcheck == F){
    stop("link must be either 'log' or 'identity'")
  }


  seasonalitycheck <- c('standard','schoolterm','none')
  if(seasonality %in% seasonalitycheck == F){
    stop("seasonality must be either 'standard' or 'schoolterm' or 'none'")
  }


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
        print('gaussian regressian failed -- switching to loess regression')
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

  if(seasonality == 'none'){

    period <- rep(1,nrow(data)-1)
    period[nrow(data)-1] <- 2

  }

  Inew <- tail(Iadjusted,-1)+1
  lIminus <- log(head(Iadjusted,-1)+1)
  Zminus <- head(Z,-1)

  pop <- data$pop

  minSmean <- max(0.01*pop,-(min(Z) - 1))
  Smean <- seq(minSmean, 0.4*mean(pop), length=250)

  loglik <- rep(NA, length(Smean))

 if(link == 'identity'){
   Inew <- log(Inew)
 }
 Inew <- round(Inew)

 if(length(input.alpha) == 0 && length(input.sbar) == 0){

   for(i in 1:length(Smean)){
     lSminus <- log(Smean[i] + Zminus)

     glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIminus) + offset(lSminus),
                   family=eval(parse(text=family))(link=link))

     loglik[i] <- glmfit$deviance

   }

   sbar <- Smean[which.min(loglik)]

   lSminus <- log(sbar + Zminus)

   glmfit <- glm(Inew ~ -1 +as.factor(period)+ (lIminus) + offset(lSminus),
                 family=eval(parse(text=family))(link=link))


   beta <- exp(head(coef(glmfit),-1))
   alpha <- tail(coef(glmfit),1)
 }


 if(length(input.alpha) == 1 && length(input.sbar) == 0){


   for(i in 1:length(Smean)){
     lSminus <- log(Smean[i] + Zminus)

     glmfit <- glm(Inew ~ -1 +as.factor(period) + offset(alpha*lIminus) + offset(lSminus),
                   family=eval(parse(text=family))(link=link))


     loglik[i] <- glmfit$deviance

   }

   sbar <- Smean[which.min(loglik)]

   lSminus <- log(sbar + Zminus)

   glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIminus) + offset(lSminus),
                 family=eval(parse(text=family))(link=link))


   beta <- exp(coef(glmfit))
 }


 if(length(input.alpha) == 0 && length(input.sbar) == 1){

   sbar <- sbar * mean(pop)
   lSminus <- log(sbar + Zminus)


   glmfit <- glm(Inew ~ -1 +as.factor(period) + (lIminus) + offset(lSminus),
                 family=eval(parse(text=family))(link=link))


   beta <- exp(head(coef(glmfit),-1))
   alpha <- tail(coef(glmfit),1)
 }


 if(length(input.alpha) == 1 && length(input.sbar) == 1){

   sbar <- sbar * mean(pop)
   lSminus <- log(sbar + Zminus)

   glmfit <- glm(Inew ~ -1 +as.factor(period)+ offset(alpha*lIminus) + offset(lSminus),
                 family=eval(parse(text=family))(link=link))

   beta <- exp(coef(glmfit))

 }

  if(seasonality == 'none'){
    beta[2] <- beta[1]
    beta <- mean(beta)
    period <- rep(1,nrow(data)-1)
  }

 confinterval <- suppressMessages(confint(glmfit))
 continterval <- confinterval[1:length(unique(period)),]
 betalow <- exp(confinterval[,1])
 betahigh <- exp(confinterval[,2])

 glmAIC <- AIC(glmfit)

 contact <- as.data.frame(cbind('time'=seq(1,length(beta[period]),1),
                                betalow[period],beta[period],betahigh[period]),row.names=F)
 names(contact) <- c('time','betalow','beta','betahigh')
 contact <- head(contact,52/IP)

  return(list('X'=X,'Y'=Y,'Yhat'=Yhat,'Smean'=Smean,'contact'=contact,'period'=period,'IP'=IP,
              'beta'=beta,'rho'=adj.rho,'Z'=Z,'pop'=pop,'time'=data$time,'AIC'=glmAIC,
              'sbar'=sbar,'alpha'=alpha,'loglik'=loglik))

}
