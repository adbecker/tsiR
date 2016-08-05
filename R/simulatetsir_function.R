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

  
  nsample <- 20
  
  inits.grid <- expand.grid(
    S0 = seq(0.01*mean(pop), 0.2*mean(pop), length=nsample), 
    I0 = seq(0.01*1e-3*mean(pop), 1*1e-3*mean(pop), length=nsample)
  )
  
  inits.res <- rep(NA,nsample*nsample)
  
  S <- rep(0,length(data$cases))
  I <- rep(0,length(data$cases))
  
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
        t0s <- epitimes(data,threshold)
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
    
  res <- matrix(0,length(data$cases),nsim)
  for(ct in 1:nsim){

    S <- rep(0,length(data$cases))
    I <- rep(0,length(data$cases))
    S[1] <- inits[[1]]
    I[1] <- inits[[2]]

    for (t in 2:(nrow(data))){

      if(pred == 'step-ahead'){
        I <- (adj.rho*data$cases)^alpha
      }
      if(pred == 'forward'){
        I <- I
      }

      lambda <- min(S[t-1],unname(beta[period[t-1]] * S[t-1] * (I[t-1])^alpha))
      
      if(is.nan(lambda) == T){lambda <- 0}
    
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
      S[t] <- max(S[t-1] + data$births[t-1] - I[t] + rnorm(n=1,mean=0,sd=add.noise.sd),0)
      #S[t] <- Z[t]+sbar + rnorm(n = 1, mean = 0, sd=add.noise.sd)
      if(S[t] < 0 && (S[t-1] + data$births[t-1] - I[t]) >0){
        warning('susceptible overflow  -- reduce additive noise sd')
      }
    }
    res[,ct] <- I / adj.rho

  }

  res[is.nan(res)] <- 0
  res[res < 1] <- 0

  res <- as.data.frame(res)
  res$mean <- rowMeans(res,na.rm=T)
  res$sd   <- apply(res, 1, function(row) sd(row[-1],na.rm=T))
  res$time <- data$time
  res$cases <- data$cases

  return(list(
    'X'=X,'Y'=Y,'Yhat' =Yhat,'pop'=pop,'Smean'=Smean,
    'beta'=beta,'rho'=adj.rho,
    'Z'=Z,'sbar'=sbar,'alpha'=alpha,'pop'=pop,
    'alphalow'=alphalow,'alphahigh'=alphahigh,
    'res'=res,'loglik'=loglik,
    'nsim'=nsim,
    'inits.grid'=inits.grid,'inits'=inits))

}
