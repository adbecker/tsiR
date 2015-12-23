#' derivative Function
#' 
#' This function computes an 8 point derivative
#' Differentiate Y with respect to X

derivative <- function(X,Y){
  weights <- c(1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280)
  
  x <- seq(X[1], X[(length(X))], length=1e4)
  y <- approxfun(X, Y)(x)
  y[1] <- 2*y[2] - y[3]
  
  dx <- diff(x)[1]
  dy <- rep(0,length(y)-8)
  
  for(i in 5:(length(y)-4)){
    dy[i-4] <- sum(weights * y[(i-4):(i+4)]) / dx  
  }
  
  dy <- approxfun(linspace(X[1], X[length(X)], length(dy)), dy)(X)
  #dy <- approx(dy,n=length(X))$y
  return(dy)
}

