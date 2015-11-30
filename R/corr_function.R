#' corr function
#' 
#' @param sim, the dataframe produced by the 'runtsir' function

corr <- function(sim){
  
  
  if(class(sim) == "list"){  
    sim <- sim$res
  }
  if(class(sim) == "data.frame"){
    sim <- sim
  }
  
  obs <- sim$cases
  pred <- sim$mean
  fit <- lm(pred ~ obs)
  c1 <- ggreg(fit)+xlab('observed')+ylab('predicted')
  print(c1)
}

