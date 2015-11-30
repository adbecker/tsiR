plotcases <- function(data){  
  
  p <- ggplot(data=data, aes(time))+
    geom_line(aes(y=cases),colour='dodgerblue',size=1) + 
    theme_bw() + theme(legend.position = "none") +
    xlab('year')+ylab('cases')
  
  
  return(p)
}

