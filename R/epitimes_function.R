#' @title epitimes
#'
#' @description The times at which we declare a new outbreak has started based on the threshold parameter.
#'
#' @param data The inputed data frame with the cases vector. This is the same data you put into runtsir.
#' @param threshold The required number of cases observed to declare it an outbreak.

epitimes <- function(data,threshold){

  x <- data$cases
  xc <- x
  xc[xc < threshold] <- 0
  index <- which(xc>0)
  test <- diff(c(0,index),lag=1)
  epitimes <- index[which(test > 1)]
  return(epitimes)
}
