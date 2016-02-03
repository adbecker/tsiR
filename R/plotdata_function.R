#' plotdata Function
#' plots the cases data as well as birth and population dynamics
#' @param data is the dataframe 

plotdata <- function(data){

  p <- ggplot(melt(data,id='time'),aes_string(x='time',y='value'))+geom_line(colour='dodgerblue',size=1)+
    facet_wrap(~variable,ncol=1,scales="free_y")+
    theme_bw()

  return(p)
}
