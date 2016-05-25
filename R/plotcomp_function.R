#' @title plotcomp
#' @description Plots just the comparison of the forward simulation fit to the data.
#' @param sim is list produced by runtsir or mcmctsir
#' @param errtype is the type of error bands to show. Defaults to '95' for 95 percent CI, the other option is 'sd' to standard deviation.
#' @param max.plot the number of individual stochastic simulations to plot. Defaults to 10. 


plotcomp <- function(sim,errtype='95',max.plot=10){

  nsim <- sim$nsim
  
  if(class(sim) == "list"){
    sim <- sim$res
  }
  if(class(sim) == "data.frame"){
    sim <- sim
  }
  
  
  drops <- c('mean','sd','error','cases','time')
  
  sim.only <- sim[,!(names(sim) %in% drops)]
  
  n <- ncol(sim.only)
  error <- qt(0.975,df=n-1)*sim$sd/sqrt(n)
  sim$error <- error

  if(errtype == '95'){
    eb <- aes(ymax = mean +  error, ymin = mean -  error)
  }
  if(errtype == 'sd'){
    eb <- aes(ymax = mean + sd, ymin = mean - sd)
  }

  comp1 <- ggplot(data=sim, aes_string('time')) + theme(legend.position = "none") +
    geom_line(aes_string(y = 'cases'), colour = "dodgerblue",size=1) + xlab('year')+ylab('cases')+
    geom_line(aes_string(y = 'mean'), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
    theme_bw()

  inversecases <- sim
  inversecases$cases <- -sim$cases

  comp2 <- ggplot(data=inversecases, aes_string('time')) + theme(legend.position = "none") +
    geom_line(aes_string(y = 'cases'), colour = "dodgerblue",size=1) + xlab('time')+ylab('cases')+
    geom_line(aes_string(y = 'mean'), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
    theme_bw()

 
  if(nsim > max.plot){
    sampledat<- sample(sim.only,max.plot) 
    sampledat$time <- sim$time
  }else{
    sampledat <- sim.only
    sampledat$time <- sim$time
  }
  
  meltdf <- melt(sampledat,id='time')
  
  comp3 <- ggplot(meltdf,aes_string(x='time',y='value',fill='variable'))+
    geom_line(alpha=0.6,colour='orangered4')+xlab('time')+ylab('cases')+
    geom_line(data=sim,aes_string(x='time',y='cases',fill=NA),colour='dodgerblue',size=1)+
    theme_bw()


  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(comp3, vp = vplayout(1, 1))
  print(comp2, vp = vplayout(2, 1))


}
