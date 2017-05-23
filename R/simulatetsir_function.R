#' simulatetsir
#' @description This function just simulates the forward prediction given the data and a parms list generated from estpars or mcmcestpars.
#'
#' @param data The data frame containing cases and interpolated births and populations.
#' @param nsim The number of simulations to do. Defaults to 100.
#' @param IP The infectious period. Defaults to 2.
#' @param parms Either the parameters estimated by estpars or mcmcestpars, or a list containing
#' beta, rho, Z, sbar, alpha, X, Y, Yhat, contact, alphalow, alphahigh, loglik, pop vectors.
#' @param method The type of next step prediction used. Options are 'negbin' for negative binomial,
#' 'pois' for poisson distribution, and 'deterministic'. Defaults to 'deterministic'.
#' @param epidemics The type of data splitting. Options are 'cont' which doesn't split the data up at all,
#' and 'break' which breaks the epidemics up if there are a lot of zeros. Defaults to 'cont'.
#' @param pred The type of prediction used. Options are 'forward' and 'step-ahead'. Defaults to 'forward'.
#' @param threshold The cut off for a new epidemic if epidemics = 'break'. Defaults to 1.
#' @param add.noise.sd The sd for additive noise, defaults to zero.
#' @param mul.noise.sd The sd for multiplicative noise, defaults to zero.
#' @param inits.fit Whether or not to fit initial conditions using simple least squares as well. Defaults to FALSE. This parameter is more necessary in more chaotic locations.


simulatetsir <- function(data, nsim = 100, IP=2,
                         parms,
                         method='deterministic',
                         epidemics='cont', pred ='forward',
                         threshold=1,
                         inits.fit=FALSE,
                         add.noise.sd = 0, mul.noise.sd = 0
){



  nzeros <- length(which(data$cases==0))
  ltot <- length(data$cases)
  if(nzeros > 0.3 * ltot && epidemics == 'cont'){
    print(sprintf('time series is %.0f%% zeros, consider using break method',100*nzeros/ltot))
  }
  Smean <- parms$Smean
  beta <- parms$beta
  adj.rho <- parms$rho
  Z <- parms$Z
  sbar <- parms$sbar
  alpha <- parms$alpha
  X <- parms$X
  Y <- parms$Y
  Yhat <- parms$Yhat
  contact <- parms$contact
  alphalow <- parms$alphalow
  alphahigh <- parms$alphahigh
  loglik <- parms$loglik
  pop <- data$pop

  datacopy <- data

  period <- rep(1:(52/IP), round(nrow(data)+1))[1:(nrow(data)-1)]

  if(IP == 1){

    period <- rep(1:(52/2),each=2, round(nrow(data)+1))[1:(nrow(data)-1)]

  }

  
  
  S <- rep(0,length(data$cases))
  I <- rep(0,length(data$cases))
  
  nsample <- 30
  
  inits.grid <- expand.grid(
    S0 = seq(0.01*mean(pop), 0.1*mean(pop), length=nsample), 
    I0 = seq(0.01*1e-3*mean(pop), 1*1e-3*mean(pop), length=nsample)
  )
  
  if(inits.fit == TRUE){
    
    inits.res <- rep(NA,nsample*nsample)
    
    for(it in 1:nrow(inits.grid)){
      S0 <- inits.grid[it,1]
      I0 <- inits.grid[it,2]
      
      S[1] <- S0
      I[1] <- I0
      
      for (t in 2:(nrow(data))){  
        
        lambda <- min(S[t-1],unname(beta[period[t-1]] * S[t-1] * (I[t-1])^alpha))
        
        #if(lambda < 1 || is.nan(lambda) == T){lambda <- 0}
        if(is.nan(lambda) == T){lambda <- 0}
        
        I[t] <- lambda 
        
        if(epidemics == 'cont'){
          I[t] <- I[t]
        }
        if(epidemics == 'break'){
          t0s <- epitimes(data,threshold)$start
          if(t %in% t0s){
            I[t] <- adj.rho[t]*data$cases[t]
          }
        }
        S[t] <- max(S[t-1] + data$births[t-1] - I[t],0)
      }
      
      inits.res[it] <- sum((I - data$cases*adj.rho)^2)
      
    }
    
    inits <- inits.grid[which.min(inits.res),]
    
    inits.grid$S0 <- inits.grid$S0/mean(pop)
    inits.grid$I0 <- inits.grid$I0/mean(pop)
    inits.grid$log10LS <- log10(inits.res)
    
    S_start <- inits[[1]]
    I_start <- inits[[2]]
    
  }else{
    
    S_start <- sbar + Z[1]
    I_start <- adj.rho[1]*datacopy$cases[1]
    
  }
  
  IC <- c(S_start,I_start)
  
  print(c('alpha'=unname(signif(alpha,2)),
          'mean beta'=unname(signif(mean(beta),3)),
          'mean rho' =unname(signif(mean(1/adj.rho),3)),
          'mean sus' =unname(signif(sbar,3)),
          'prop. init. sus.' =unname(signif(S_start/mean(pop),3)),
          'prop. init. inf.' =unname(signif(I_start/mean(pop),3))    
  ))
  
  
  nsim <- nsim
  res <- matrix(0,length(data$cases),nsim)
  for(ct in 1:nsim){
    
    S <- rep(0,length(data$cases))
    I <- rep(0,length(data$cases))
    
    S[1] <- S_start
    I[1] <- I_start
    
    for (t in 2:(nrow(data))){
      
      
      if(pred == 'forward'){
        I <- I
      }
      
      lambda <- min(S[t-1],unname(beta[period[t-1]] * S[t-1] * (I[t-1])^alpha))
      
      if(pred == 'step-ahead'){
        lambda <- min(S[t-1],unname(beta[period[t-1]] * S[t-1] * (adj.rho[t]*data$cases[t])^alpha))
      }
      
      #if(lambda < 1 || is.nan(lambda) == T){lambda <- 0}
      if(is.nan(lambda) == T){lambda <- 0}
      
      if(method == 'deterministic'){
        I[t] <- lambda * rnorm( n = 1, mean = 1, sd=mul.noise.sd)
        if(I[t] < 0 && lambda >= 0 ){
          warning('infected overflow  -- reduce multiplicative noise sd')
        }
      }
      if(method == 'negbin'){
        I[t] <- rnbinom(n=1,mu=lambda,size=I[t-1]+1e-10)
      }
      if(method == 'pois'){
        I[t] <- rpois(n=1,lambda=lambda)
      }
      if(epidemics == 'cont'){
        I[t] <- I[t]
      }
      if(epidemics == 'break'){
        
        t0s <- epitimes(data,threshold)$start
        if(t %in% t0s){
          I[t] <- adj.rho[t]*data$cases[t]
        }
      }
      S[t] <- max(S[t-1] + data$births[t-1] - I[t] + rnorm(n=1,mean=0,sd=add.noise.sd),0)
      
      if(S[t] < 0 && (S[t-1] + data$births[t-1] - I[t]) >0 ){
        warning('susceptible overflow  -- reduce additive noise sd')
      }
    }
    res[,ct] <- I / adj.rho
    
  }
  
  res[is.nan(res)] <- 0
  res[res < 1] <- 0
  
  res <- as.data.frame(res)
  
  #res$mean <- apply(res, 1, function(row) mean(row[-1],na.rm=T))
  res$mean <- rowMeans(res,na.rm=T)
  res$sd   <- apply(res, 1, function(row) sd(row[-1],na.rm=T))
  res$time <- data$time
  res$cases <- data$cases
  
  obs <- res$cases
  pred <- res$mean
  
  fit <- lm(pred ~ obs)
  
  rsquared <- signif(summary(fit)$adj.r.squared, 2)
  
  

  return(list(
    'X'=X,'Y'=Y,'Yhat' =Yhat,'pop'=pop,'Smean'=Smean,
    'beta'=beta,'rho'=adj.rho,
    'Z'=Z,'sbar'=sbar,'alpha'=alpha,'pop'=pop,
    'alphalow'=alphalow,'alphahigh'=alphahigh,
    'res'=res,'loglik'=loglik,
    'nsim'=nsim,'contact'=contact,'rsquared'=rsquared,
    'inits.fit'=inits.fit,
    'inits.grid'=inits.grid,'inits'=IC))

}
