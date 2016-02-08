#' @title residuals.births
#' @description computes the residuals for when X is the cumulative births. Used internally.
#' @param rho The reporting rate, used to get units correct.
#' @param Yhat The fitted regression line.
#' @param Y The cumulative cases.

residual.births <- function(rho,Yhat,Y){
  ## when X is cumbirths
  Z <- -rho*(Y-Yhat)
  return(Z)
}
