#' @title plotcases
#' @description Plots just the cases data.
#' @param data The data frame with cases.
plotcases <- function(data){

  p <- ggplot(data=data, aes_string('time'))+
    geom_line(aes_string(y='cases'),colour='dodgerblue',size=1) +
    theme_bw() + theme(legend.position = "none") +
    xlab('year')+ylab('cases')


  return(p)
}
