jagsfilter <- function(mcmcresults){
  
  probabilities <- c(0.025, 0.975)
  
  pull <- t(apply(do.call(rbind, mcmcresults), 2, function(f) {
    c(mean=mean(f), sd=sd(f), quantile(f, probabilities))
  }))
  
  
  pull <- pull[,c('mean', '2.5%', '97.5%')]
  return(pull)
}



