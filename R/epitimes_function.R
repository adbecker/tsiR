#' epitimes Function
#'
#' @param data dataframe with the cases vector.
#' @param threshold the required number of cases to call an outbreak.

epitimes <- function(data,threshold){

  x <- data$cases
  xc <- x
  xc[xc < threshold] <- 0
  index <- which(xc>0)
  test <- diff(c(0,index),lag=1)
  epitimes <- index[which(test > 1)]
  return(epitimes)
}
