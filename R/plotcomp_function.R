
plotcomp <- function(sim,errtype='95'){
  
  if(class(sim) == "list"){  
    sim <- sim$res
  }
  if(class(sim) == "data.frame"){
    sim <- sim
  }
  n <- nrow(sim)
  error <- qt(0.975,df=n-1)*sim$sd/sqrt(n)
  sim$error <- error

  if(errtype == '95'){
    eb <- aes(ymax = mean +  error, ymin = mean -  error)
  }
  if(errtype == 'sd'){
    eb <- aes(ymax = mean + sd, ymin = mean - sd)
  }

  comp1 <- ggplot(data=sim, aes(time)) + theme(legend.position = "none") +
    geom_line(aes(y = cases), colour = "dodgerblue",size=1) + xlab('year')+ylab('cases')+
    geom_line(aes(y = mean), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
    theme_bw()

  comp2 <- ggplot(data=sim, aes(time)) + theme(legend.position = "none") +
    geom_line(aes(y = -cases), colour = "dodgerblue",size=1) + 
    geom_line(aes(y = mean), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
    theme_bw()+xlab('time')+ylab('cases')

  meltdf <- melt(subset(sim,select=-c(mean,sd,error)),id='time')
  
  comp3 <- ggplot(meltdf,aes(x=time,y=value))+
    geom_line(alpha=0.6,colour='orangered4')+xlab('time')+ylab('cases')+
    geom_line(data=sim,aes(x=time,y=cases),colour='dodgerblue',size=1)+
    theme_bw()
  
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(comp3, vp = vplayout(1, 1))
  print(comp2, vp = vplayout(2, 1))


}

