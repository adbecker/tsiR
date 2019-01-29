#' @title residuals.cases
#' @description Computes the residuals for when X is the cumulative cases. Used internally.
#' @param Yhat The fitted regression line.
#' @param Y The cumulative births.


residual.cases <- function(Yhat,Y){
  ## when X is cumcases can do standard calculation
  Z <- Y - Yhat
  return(Z)
}
