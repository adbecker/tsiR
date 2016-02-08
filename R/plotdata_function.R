#' @title plotdata
#' @description Plots the cases data as well as birth and population dynamics.
#' @param data The dataframe with time, cases, births, and pop.

plotdata <- function(data){

  p <- ggplot(melt(data,id='time'),aes_string(x='time',y='value'))+geom_line(colour='dodgerblue',size=1)+
    facet_wrap(~variable,ncol=1,scales="free_y")+
    theme_bw()

  return(p)
}
