#' plotres function
#' 
#' function to plot diagnostics and results of the runtsir function.
#' @param dat, the list produced from the runtsir function or the mcmctsir function
#' 
plotres <- function(dat){
  regdf <- NULL
  regdf$X <- dat$X
  regdf$Y <- dat$Y
  regdf$Yhat <- dat$Yhat
  regdf$time <- dat$res$time
  regdf <- as.data.frame(regdf)
  meltdf <- melt(regdf, id.vars = "time")
  
  p1 <- ggplot(meltdf, aes(time, value,col=variable)) + geom_line(size=2)+theme_bw()
  
  rhodf <- NULL
  rhodf$time <- dat$res$time
  rhodf$rho <- dat$rho
  rhodf <- as.data.frame(rhodf)
  
  p2 <- ggplot(rhodf,aes(time,1/rho))+geom_line(size=2)+theme_bw()+
    ggtitle(sprintf('mean rho=%g',signif(mean(1/dat$rho),2)))
  
  
  
  resdf <- NULL
  resdf$time <- dat$res$time
  resdf$Z <- dat$Z
  resdf$S <- dat$Z + dat$sbar
  resdf <- as.data.frame(resdf)
  meltres <- melt(resdf, id.vars = "time")
  
  p3 <- ggplot(meltres, aes(time, value,col=variable)) + geom_line(size=2)+theme_bw()
  
  loglikdf <- NULL
  loglikdf$sbar <- dat$Smean
  loglikdf$loglik <- dat$loglik
  loglikdf <- as.data.frame(loglikdf)
  p9 <- ggplot(loglikdf,aes(sbar,loglik))+geom_line()+geom_point()+
    theme_bw()+geom_vline(xintercept = dat$sbar,linetype = "longdash")+
    ggtitle(sprintf('Sbar=%g: %g of population',signif(dat$sbar,2),signif(dat$sbar/mean(dat$pop),2)))
  
  betadf <- NULL
  betadf$time <- seq(1,length(dat$beta),1)
  betadf$beta <- dat$beta
  betadf <- as.data.frame(betadf)
  
  
  p4 <- ggplot(betadf,aes(time,beta))+geom_line(size=2)+theme_bw()+
    ggtitle(sprintf('mean beta = %g, alpha=%g',signif(mean(dat$beta),2),signif(dat$alpha,2)))
  
  
  if('contact' %in% names(dat)){
    
    p4 <- ggplot(dat$contact,aes(time,beta))+geom_line(size=2)+
      geom_ribbon(ymin=dat$contact$betalow,ymax=dat$contact$betahigh,alpha=0.5,col='dodgerblue',fill='dodgerblue')+
      ylim(c(min(dat$contact$betalow),max(dat$contact$betahigh)))+theme_bw()+
      ggtitle(sprintf('mean beta = %g, alpha=%g',signif(mean(dat$beta),2),signif(dat$alpha,2)))
    
    
  }
  
  obs <- dat$res$cases
  pred <- dat$res$mean
  fit <- lm(pred ~ obs)
  
  p5 <- ggreg(fit)+xlab('observed')+ylab('predicted')
  
  
  
  n <- nrow(dat$res)
  error <- qt(0.975,df=n-1)*dat$res$sd/sqrt(n)
  dat$res$error <- error
  
  eb <- aes(ymax = mean +  error, ymin = mean -  error)
  
  p6 <- ggplot(data=dat$res, aes(time)) + theme(legend.position = "none") +
    geom_line(aes(y = cases), colour = "dodgerblue",size=1) + xlab('year')+ylab('cases')+
    geom_line(aes(y = mean), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
    theme_bw()
  
  p7 <- ggplot(data=dat$res, aes(time)) + theme(legend.position = "none") +
    geom_line(aes(y = -cases), colour = "dodgerblue",size=1) + xlab('time')+ylab('cases')+
    geom_line(aes(y = mean), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
    theme_bw()
  
  meltdf <- melt(subset(dat$res,select=-c(mean,sd,error)),id='time')
  
  p8 <- ggplot(meltdf,aes(x=time,y=value))+
    geom_line(alpha=0.6,colour='orangered4')+ylab('cases')+xlab('time')+
    geom_line(data=dat$res,aes(x=time,y=cases),colour='dodgerblue',size=1)+
    theme_bw()
  
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, 2)))
  
  if(all(is.na(dat$loglik)) == T){
  print(p1,vp=vplayout(1,1))
  print(p2,vp=vplayout(1,2))
  print(p3,vp=vplayout(2,1:2))
  print(p4,vp=vplayout(3,1))
  print(p5,vp=vplayout(3,2))
  print(p8,vp=vplayout(4,1:2))
  print(p7,vp=vplayout(5,1:2))
  }
  
  
  print(p1,vp=vplayout(1,1))
  print(p2,vp=vplayout(1,2))
  print(p3,vp=vplayout(2,1))
  print(p9,vp=vplayout(2,2))
  print(p4,vp=vplayout(3,1))
  print(p5,vp=vplayout(3,2))
  print(p8,vp=vplayout(4,1:2))
  print(p7,vp=vplayout(5,1:2))
  
  
  
  
}