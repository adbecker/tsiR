residual.cases <- function(Yhat,Y){
  ## when X is cum.cases
  Z <- Y - Yhat
  return(Z)
}

