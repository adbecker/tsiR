tsiRdata <- function(time,cases,births,pop,IP=2){
  
  intcases <- as.numeric(tapply(cases, (seq_along(cases)-1) %/% IP, sum))
  intweek <- time[seq(1, length(time), IP)]
  intbirths <- approx(births,n=length(intweek))$y / (52 / IP)
  intpop <- approx(pop,n=length(intweek))$y 
  
  data <- as.data.frame(cbind(intweek,intcases,intbirths,intpop))
  names(data) <- c('time','cases','births','pop')
  
  return(data)
}