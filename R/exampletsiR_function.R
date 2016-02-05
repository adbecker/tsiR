#' example tsiR function

exampletsiR <- function(){
  
  cat("\nFirst load data -- here we will use the London data")
  data <- twentymeas[["London"]]
  
  cat("\nDesignate an IP (infectious period) in weeks -- here we will use 2")
  IP <- 2
  
  cat('\nHow many simulations do you want to do?')
  nsim <- scan(what=numeric(),nmax=1,quiet = T)
  
  cat('\nAssign a value to alpha (just hit enter if you want to estimate)')
  alpha <- scan(what=numeric(),nmax=1,quiet = T)
  
  cat('\nAssign a value to sbar (between 0 and 1) (just hit enter if you want to estimate)')
  sbar <- scan(what=numeric(),nmax=1,quiet = T)
  
  cat('\nDo you want births or cases on the x axis? Options are "cumcases" or "cumbirths"')
  xreg <- scan(what=character(),nmax=1,quiet = T)
  
  cat('\nWhat type of regression do you want to use? Options are 
"gaussian", "lm" (linear model), "spline" (smooth.spline with 2.5 degrees freedom),
"lowess" (with f = 2/3, iter = 1), "loess" (degree 1), and "user". Default is "gaussian". 
    "user" is a user inputted vector')
  regtype <- scan(what=character(),nmax=1,quiet = T)
  
  if(regtype == 'user'){
    
    cat('\nPlease enter the user regression function')
    userYhat <- scan(what=numeric(),nmax=1,quiet = T)
    
  }
  
  cat('\nWhat type of regression family? Options are "gaussian" and "poisson"')
  family <- scan(what=character(),nmax=1,quiet = T)
  
  cat('\nWhat type of seasonality? Options are "standard" for the 52/IP point contact function
    and "schoolterm" for on off of the school calender')
  seasonality <- scan(what=character(),nmax=1,quiet = T)
  
  cat('\nWhat type of forward simulation distribution? Options are "negbin" and "pois" and "deterministic"')
  method <- scan(what=character(),nmax=1,quiet = T)
  
  if(method == 'deterministic'){
    
    cat('\nDo you want to have multiplicative noise? This determines the standard deviation of the noise')
    mul.noise.sd <-  scan(what=numeric(),nmax=1,quiet = T)
    
    cat('\nDo you want to have additive noise? This determines the standard deviation of the noise')
    add.noise.sd <-  scan(what=numeric(),nmax=1,quiet = T)
    
  }else{
    
    mul.noise.sd <- 0
    add.noise.sd <- 0
    
  }
  
  #cat('\nWhat type of forward simulation? Options are "cont" and "break". For London we use "cont"')
  #epidemics <- scan(what=character(),nmax=1,quiet = T)
  epidemics <- 'cont'
  
  cat('\nWhat type of forward prediction? Options are "forward" and "step-ahead"')
  pred <- scan(what=character(),nmax=1,quiet = T)
  
  cat('\nRunning the code...')
  
  res <- runtsir(data=data,xreg=xreg,IP=IP,nsim=nsim,
                 regtype=regtype,userYhat=userYhat,alpha=alpha,
                 sbar=sbar,family=family,method=method,epidemics=epidemics,
                 pred=pred,seasonality=seasonality,
                 add.noise.sd=add.noise.sd,mul.noise.sd=mul.noise.sd
  )
  
  if(res$rsquared < 0.3){
    cat("\nThe fit doesn't seem that great -- try different options?")
  }
  
  cat('\nDo you want to plot? "yes" or "no"')
  ploton <- scan(what=character(),nmax=1,quiet = T)
  
  if(ploton == 'yes'){
    plotres(res)
  }
  
  return(res)
  
}
