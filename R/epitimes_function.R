#' @title epitimes
#'
#' @description The times at which we declare a new outbreak has started based on the threshold parameter.
#'
#' @param data The inputed data frame with the cases vector. This is the same data you put into runtsir.
#' @param threshold The required number of cases observed to declare it an outbreak.
#' @param epi.length The required duration (in 52/IP weeks)  to declare it an outbreak.

epitimes <- function(data,threshold,epi.length=3){
  
  index <- which(data$cases[1:nrow(data) - 1] > threshold & data$cases[2:nrow(data)] > threshold)
  index.length <- length(index)
  
  start <- ifelse(diff(index,lag=1) > 1, index[2:index.length], NA)
  start <- c(index[1], start)
  start <- start[!is.na(start)]
  
  end <- ifelse(diff(index,lag=1) > 1, index[1:index.length - 1], NA)
  end <- c(end, index[index.length])
  end <- end[!is.na(end)]
  
  index.diff <- which(start[2:length(start)] - end[1:length(end) - 1] <= epi.length)
  
  for(i in 1:length(index.diff)){
    new.i <- index.diff[i]
    start[new.i + 1] <- NA
    end[new.i] <- ifelse(is.na(end[new.i]) == FALSE, end[new.i + 1], NA)
    end[new.i + 1] <- NA
  }
  
  start <- start[!is.na(start)]
  end <- end[!is.na(end)]
  
  duration <- end - start
  
  final.size <- rep(NA,length(start))
  peak.size <- rep(NA,length(start))
  if(length(start) > 0){
    for (i in 1:length(start)) {
      final.size[i] <- sum(data$cases[start[i]:end[i]])
      peak.size[i] <- max(data$cases[start[i]:end[i]])
    }
  }
  epi <- as.data.frame(cbind(start, end, duration,final.size,peak.size))
  epi <- epi[!is.na(epi$duration), ]
  epi <- epi[epi$duration > epi.length, ] #only keep epidemics that are greater than a specified length
  
  return(epi)
}


