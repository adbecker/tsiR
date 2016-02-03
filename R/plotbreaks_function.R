#' plotbreaks function
#' 
#' @param data, dataframe with the cases vector.
#' @param threshold, the epidemic threshold.
plotbreaks <- function(data,threshold){
  t0s <- epitimes(data,threshold) 
  ptsdat <- as.data.frame(cbind('t'=data$time[t0s],'cpts'=data$cases[t0s]))
  p <- ggplot()+
    geom_line(data=data, aes_string(x='time',y='cases'),colour='dodgerblue',size=1) + 
    theme_bw() + theme(legend.position = "none") +
    xlab('time')+ylab('cases') + 
    geom_vline(xintercept = data$time[t0s],linetype='dashed')
  #+geom_point(data=ptsdat,aes(x=t,y=cpts),col='orangered1',size=3)
  
  return(p) 
}

