#' residuals.births Function
#' computes the residuals for when X is the cumulative births
#' used internally
#' @param rho is the reporting rate to get units correct
#' @param Yhat is the fitted regression
#' @param Y is the cumulative cases

residual.births <- function(rho,Yhat,Y){
  ## when X is cumbirths
  Z <- -rho*(Y-Yhat)
  return(Z)
}
