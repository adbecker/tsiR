#' @title tsiRdata
#' @description A function to take in time cases births and pop vectors (of any lengths) and interpolate them using the given infectious period.
#' @param time The time vector.
#' @param cases The cases vector.
#' @param births The births vector.
#' @param pop The population vector.
#' @param IP The infectious period (in weeks) to discretize to. Defaults to 2.

tsiRdata <- function(time,cases,births,pop,IP=2){

  ## sum cases by IP weeks
  intcases <- as.numeric(tapply(cases, (seq_along(cases)-1) %/% IP, sum))

  ## select only IP week weeks
  intweek <- time[seq(1, length(time), IP)]

  ## interpolate yearly births and pop by IP intervals
  intbirths <- approx(births,n=length(intweek))$y / (52 / IP)
  intpop <- approx(pop,n=length(intweek))$y

  ## combine into a data frame with correct names
  data <- as.data.frame(cbind(intweek,intcases,intbirths,intpop))
  names(data) <- c('time','cases','births','pop')

  return(data)
}
