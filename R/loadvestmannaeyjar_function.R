#' loadvestmannaeyjar function
#' 
#' function to load data from vestmannaeyjar, Iceland.
#' the period here is 24 weeks, therefore IP should be (52/24)

loadvestmannaeyjar <- function(){
  
  load('vestmanhaeyjar.RData') 
  names(data) <- c('time','births','pop','cases')
  data <- subset(data,time < 1965)
  return(data)
  
}