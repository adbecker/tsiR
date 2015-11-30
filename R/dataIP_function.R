dataIP <- function(data,demog,IP){
  
  bicases <- as.numeric(tapply(data$cases, (seq_along(data$cases)-1) %/% IP, sum))
  biweek <- data$time[seq(1, length(data$time), IP)]
  intbirths <- approx(demog$births,n=length(biweek))$y / (52 / IP)
  intpop <- approx(demog$pop,n=length(biweek))$y 
  
  data <- as.data.frame(cbind(biweek,bicases,intbirths,intpop))
  names(data) <- c('time','cases','births','pop')
  
  return(data)
}

