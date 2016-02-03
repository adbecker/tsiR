#' residuals.cases Function
#' computes the residuals for when X is the cumulative cases
#' used internally
#' @param Yhat is the fitted regression
#' @param Y is the cumulative births


residual.cases <- function(Yhat,Y){
  ## when X is cumcases
  Z <- Y - Yhat
  return(Z)
}
