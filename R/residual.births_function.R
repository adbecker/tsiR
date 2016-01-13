residual.births <- function(rho,Yhat,Y){
  ## when X is cum.births
  Z <- -rho*(Y-Yhat)
  return(Z)
}


