residual.births <- function(Yhat,Y){
  ## when X is cum.births
  Z <- -(Y-Yhat)
  return(Z)
}


