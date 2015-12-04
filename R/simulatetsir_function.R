simulatetsir <- function(data, nsim = 100, IP=2,
                         parms,
                         method='deterministic',
                         epidemics='cont', pred ='forward',
                         threshold=1,
                         add.noise.sd = 0, mul.noise.sd = 0
){
  
  
  
  nzeros <- length(which(data$cases==0))
  ltot <- length(data$cases)
  if(nzeros > 0.3 * ltot && epidemics == 'cont'){
    print(sprintf('time series is %.0f%% zeros, consider using break method',100*nzeros/ltot))
  }
  
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
  pop <- parms$pop
  
  datacopy <- data
  
  period <- rep(1:(52/IP), round(nrow(data)+1))[1:(nrow(data)-1)]
  
  if(IP == 1){
    
    period <- rep(1:(52/2),each=2, round(nrow(data)+1))[1:(nrow(data)-1)]
    
  }  
  
  res <- matrix(0,length(data$cases),nsim)
  for(ct in 1:nsim){
    
    S <- rep(0,length(data$cases))
    I <- rep(0,length(data$cases))
    S[1] <- sbar+Z[1]
    I[1] <- datacopy$cases[1] * adj.rho[1]
    
    for (t in 2:(nrow(data))){
      
      if(pred == 'step-ahead'){
        I <- (adj.rho*data$cases)^alpha
      }
      if(pred == 'forward'){
        I <- I
      }
      
      lambda <- unname(beta[period[t-1]] * S[t-1] * (I[t-1])^alpha)
      
      if(lambda < 1 || is.nan(lambda) == T){lambda <- 0}
      
      if(method == 'deterministic'){
        I[t] <- lambda * rnorm( n = 1, mean = 1, sd=mul.noise.sd)
        if(I[t] < 0 && lambda >= 0 ){
          warning('infected overflow  -- reduce multiplicative noise sd')
        }
      }
      if(method == 'negbin'){
        I[t] <- rnbinom(n=1,mu=lambda,size=I[t-1])
      }
      if(method == 'pois'){
        I[t] <- rpois(n=1,lambda=lambda)
      }
      if(epidemics == 'cont'){
        I[t] <- I[t]
      }
      if(epidemics == 'break'){
        
        t0s <- epitimes(data,threshold)
        if(t %in% t0s){
          I[t] <- adj.rho[t]*data$cases[t]
        }
      }
      S[t] <- Z[t]+sbar + rnorm(n = 1, mean = 0, sd=add.noise.sd)
      if(S[t] < 0 && (Z[t] + sbar) >0){
        warning('susceptible overflow  -- reduce additive noise sd')
      }
    }
    res[,ct] <- I / adj.rho
    
  }
  
  res[is.nan(res)] <- 0
  res[res < 1] <- 0
  
  res <- as.data.frame(res)
  res$mean <- apply(res, 1, function(row) mean(row[-1],na.rm=T))
  res$sd   <- apply(res, 1, function(row) sd(row[-1],na.rm=T))
  res$time <- data$time
  res$cases <- data$cases
  
  return(list(
  'X'=X,'Y'=Y,'Yhat' =Yhat,'pop'=pop,
  'beta'=beta,'contact'=contact,'rho'=adj.rho,
  'Z'=Z,'sbar'=sbar,'alpha'=alpha,
  'alphalow'=alphalow,'alphahigh'=alphahigh,
  'res'=res,'loglik'=loglik,
  'nsim'=nsim))
  
}