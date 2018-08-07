#' @title plotres
#'
#' @description Plots diagnostics and results of the runtsir function.
#' @param dat the list produced from the runtsir, mcmctsir, and simulatetsir function.
#' @param max.plot the number of individual stochastic simulations to plot. Defaults to 10.
plotres <- function(dat,max.plot = 10){

  ##note I have to use aes_string and then characters instead of aes to get around the R CMD check

  if(is.null(dat$SIRS) == TRUE){
    dat$SIRS <- FALSE
  }

  if(dat$SIRS == TRUE){

    ll.melt <- dat$ll.melt
    p1 <- ggplot(ll.melt, aes_string(x='X1',y='X2', z='loglik')) +
      geom_tile(aes_string(fill= 'loglik')) + scale_fill_gradient(low="white", high="red") + theme_bw()+
      geom_contour(col='black')+
      geom_point(aes(x=dat$k,y=dat$m),col='black')+
      xlab('strength of immunity')+
      ylab('duration of immunity')+
      ggtitle(sprintf('duration = %g, strength = %g',round(dat$m),signif(dat$k,2)))


    update.ll <- ll.melt[,which(!apply(ll.melt==0,2,all))]

    if(nrow(update.ll)*ncol(update.ll) != nrow(ll.melt)*ncol(ll.melt)){

      lldf <- NULL
      lldf$loglik <- update.ll$loglik
      lldf$Var <- update.ll[names(update.ll)[which(names(update.ll) != 'loglik')]]
      lldf <- as.data.frame(lldf)
      names(lldf) <- c('loglik','Var')

      p1 <- ggplot(lldf,aes_string('Var','loglik'))+geom_line(size=2)+theme_bw()+
        ggtitle(sprintf('duration = %g, strength = %g',round(dat$m),signif(dat$k,2)))

    }

    Sdf <- NULL
    Sdf$time <- head(dat$res$time,-1)
    Sdf$S <- dat$S
    Sdf <- as.data.frame(Sdf)

    p2 <- ggplot(Sdf,aes_string('time','S'))+geom_line(size=2)+theme_bw()+
      ggtitle(bquote(bar(rho)==.(signif(mean(1/dat$rho),2))))

    p3 <- ggplot(dat$contact,aes_string('time','beta'))+geom_line(size=2)+
      geom_ribbon(ymin=dat$contact$betalow,ymax=dat$contact$betahigh,alpha=0.5,col='dodgerblue',fill='dodgerblue')+
      ylim(c(min(dat$contact$betalow),max(dat$contact$betahigh)))+theme_bw()+
      ggtitle(bquote(bar(beta)==.(signif(mean(dat$beta),2))~','~alpha==.(signif(dat$alpha,3))))+
      ylab(bquote(beta))

    if(sum(sum(is.na(dat$contact))) > 0){

      p3 <- ggplot(betadf,aes_string('time','beta'))+geom_line(size=2)+theme_bw()+
        ggtitle(bquote(bar(beta)==.(signif(mean(dat$beta),2))~','~alpha==.(signif(dat$alpha,3))))+
        ylab(bquote(beta))


    }

    p4 <- logcorr(dat)+geom_abline(slope = 1,colour='dodgerblue')


    drops <- c('mean','sd','error','cases','time')

    sim.only <- dat$res[,!(names(dat$res) %in% drops)]

    n <- ncol(sim.only)
    error <- qt(0.975,df=n-1)*dat$res$sd/sqrt(n)
    dat$res$error <- error

    eb <- aes(ymax = mean +  error, ymin = mean -  error)

    p6 <- ggplot(data=dat$res, aes_string('time')) + theme(legend.position = "none") +
      geom_line(aes_string(y = 'cases'), colour = "dodgerblue",size=1) + xlab('year')+ylab('cases')+
      geom_line(aes_string(y = 'mean'), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
      theme_bw()

    inversecases <- dat$res
    inversecases$cases <- -dat$res$cases

    p7 <- ggplot(data=inversecases, aes_string('time')) + theme(legend.position = "none") +
      geom_line(aes_string(y = 'cases'), colour = "dodgerblue",size=1) + xlab('time')+ylab('cases')+
      geom_line(aes_string(y = 'mean'), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
      theme_bw()

    if(dat$nsim > max.plot){
      sampledat<- sample(sim.only,max.plot)
      sampledat$time <- dat$res$time
    }else{
      sampledat <- sim.only
      sampledat$time <- dat$res$time
    }

    meltdf <- melt(sampledat,id='time')

    p8 <- ggplot(meltdf,aes_string(x='time',y='value',fill='variable'))+
      geom_line(alpha=0.6,colour='orangered4')+xlab('time')+ylab('cases')+
      geom_line(data=dat$res,aes_string(x='time',y='cases',fill=NA),colour='dodgerblue',size=1)+
      theme_bw()



    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 2)))

    print(p1,vp=vplayout(1,1))
    print(p2,vp=vplayout(1,2))
    print(p3,vp=vplayout(2,1))
    print(p4,vp=vplayout(2,2))
    print(p8,vp=vplayout(3,1:2))
    print(p7,vp=vplayout(4,1:2))


  }else{



    regdf <- NULL
    regdf$X <- dat$X
    regdf$Y <- dat$Y
    regdf$Yhat <- dat$Yhat
    regdf$time <- dat$res$time
    regdf <- as.data.frame(regdf)
    meltdf <- melt(regdf, id.vars = "time")

    p1 <- ggplot(meltdf, aes_string('time', 'value',col='variable')) + geom_line(size=2)+theme_bw()

    rhodf <- NULL
    rhodf$time <- dat$res$time
    rhodf$rho <- 1/dat$rho
    rhodf <- as.data.frame(rhodf)

    p2 <- ggplot(rhodf,aes_string('time','rho'))+geom_line(size=2)+theme_bw()+
      ggtitle(bquote(bar(rho)==.(signif(mean(1/dat$rho),2))))+
      ylab(bquote(1/rho))


    resdf <- NULL
    resdf$time <- dat$res$time
    resdf$Z <- dat$Z
    resdf$S <- dat$Z + dat$sbar
    resdf <- as.data.frame(resdf)
    meltres <- melt(resdf, id.vars = "time")

    p3 <- ggplot(meltres, aes_string('time', 'value',col='variable')) + geom_line(size=2)+theme_bw()

    loglikdf <- NULL
    loglikdf$sbar <- dat$Smean
    loglikdf$loglik <- dat$loglik
    loglikdf <- as.data.frame(loglikdf)
    p9 <- ggplot(loglikdf,aes_string('sbar','loglik'))+geom_line()+geom_point()+
      theme_bw()+geom_vline(xintercept = dat$sbar,linetype = "longdash")+
      ggtitle(bquote(bar(S) ==.(signif(dat$sbar,2))~','~.(signif(dat$sbar/mean(dat$pop)*100,2))
                     ~'%'))+
      xlab(bquote(bar(S)))


    betadf <- NULL

    betadf <- dat$contact

    # betadf$time <- seq(1,length(dat$beta),1)
    # betadf$beta <- dat$beta
    betadf <- as.data.frame(betadf)


    p4 <- ggplot(betadf,aes_string('time','beta'))+geom_line(size=2)+theme_bw()+
      ggtitle(bquote(bar(beta)==.(signif(mean(dat$beta),2))~','~alpha==.(signif(dat$alpha,3))))+
      ylab(bquote(beta))


    if('contact' %in% names(dat)){

      p4 <- ggplot(dat$contact,aes_string('time','beta'))+geom_line(size=2)+
        geom_ribbon(ymin=dat$contact$betalow,ymax=dat$contact$betahigh,alpha=0.5,col='dodgerblue',fill='dodgerblue')+
        ylim(c(min(dat$contact$betalow),max(dat$contact$betahigh)))+theme_bw()+
        ggtitle(bquote(bar(beta)==.(signif(mean(dat$beta),2))~','~alpha==.(signif(dat$alpha,3))))+
        ylab(bquote(beta))

      if(sum(sum(is.na(dat$contact))) > 0){

        p4 <- ggplot(betadf,aes_string('time','beta'))+geom_line(size=2)+theme_bw()+
          ggtitle(bquote(bar(beta)==.(signif(mean(dat$beta),2))~','~alpha==.(signif(dat$alpha,3))))+
          ylab(bquote(beta))


      }

    }

    p4 <- p4 + xlab(sprintf('time mod %g',nrow(dat$contact)))

    #p5 <- logcorr(dat)+geom_abline(slope = 1,colour='dodgerblue')


    if(dat$inits.fit == TRUE){
      inits.grid <- dat$inits.grid

      p5 <-  ggplot(inits.grid, aes_string(x='S0',y='I0', z='log10LS')) +
        geom_tile(aes_string(fill= 'log10LS')) + scale_fill_gradient(low="white", high="black") + theme_bw()+
        geom_contour(col='black')+
        geom_point(aes(x=dat$inits[1]/mean(dat$pop),y=dat$inits[2]/mean(dat$pop)),col='red')+
        xlab('prop. init. sus.')+
        ylab('prop. init. inf.')
    }

    drops <- c('mean','sd','error','cases','time')

    sim.only <- dat$res[,!(names(dat$res) %in% drops)]

    n <- ncol(sim.only)
    error <- qt(0.975,df=n-1)*dat$res$sd/sqrt(n)
    dat$res$error <- error
    eb <- aes(ymax = mean +  error, ymin = mean -  error)

    p6 <- ggplot(data=dat$res, aes_string('time')) + theme(legend.position = "none") +
      geom_line(aes_string(y = 'cases'), colour = "dodgerblue",size=1) + xlab('year')+ylab('cases')+
      geom_line(aes_string(y = 'mean'), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
      theme_bw()

    inversecases <- dat$res
    inversecases$cases <- -dat$res$cases

    p7 <- ggplot(data=inversecases, aes_string('time')) + theme(legend.position = "none") +
      geom_line(aes_string(y = 'cases'), colour = "dodgerblue",size=1) + xlab('time')+ylab('cases')+
      geom_line(aes_string(y = 'mean'), colour = "orangered4",size=1) + geom_ribbon(eb,alpha=0.3)+
      theme_bw()

    if(dat$nsim > max.plot){
      sampledat<- sample(sim.only,max.plot)
      sampledat$time <- dat$res$time
    }else{
      sampledat <- sim.only
      sampledat$time <- dat$res$time
    }

    meltdf <- melt(sampledat,id='time')

    p8 <- ggplot(meltdf,aes_string(x='time',y='value',fill='variable'))+
      geom_line(alpha=0.6,colour='orangered4')+xlab('time')+ylab('cases')+
      geom_line(data=dat$res,aes_string(x='time',y='cases',fill=NA),colour='dodgerblue',size=1)+
      theme_bw()


    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5, 2)))
    if(dat$inits.fit == TRUE){
      if(all(is.na(dat$loglik)) == T){
        print(p1,vp=vplayout(1,1))
        print(p2,vp=vplayout(1,2))
        print(p3,vp=vplayout(2,1:2))
        print(p4,vp=vplayout(3,1))
        print(p5,vp=vplayout(3,2))
        print(p8,vp=vplayout(4,1:2))
        print(p7,vp=vplayout(5,1:2))
      }else{
        print(p1,vp=vplayout(1,1))
        print(p2,vp=vplayout(1,2))
        print(p3,vp=vplayout(2,1))
        print(p9,vp=vplayout(2,2))
        print(p4,vp=vplayout(3,1))
        print(p5,vp=vplayout(3,2))
        print(p8,vp=vplayout(4,1:2))
        print(p7,vp=vplayout(5,1:2))
      }
    }else{
      if(all(is.na(dat$loglik)) == T){
        print(p1,vp=vplayout(1,1))
        print(p2,vp=vplayout(1,2))
        print(p3,vp=vplayout(2,1:2))
        print(p4,vp=vplayout(3,1:2))
        print(p8,vp=vplayout(4,1:2))
        print(p7,vp=vplayout(5,1:2))
      }else{
        print(p1,vp=vplayout(1,1))
        print(p2,vp=vplayout(1,2))
        print(p3,vp=vplayout(2,1))
        print(p9,vp=vplayout(2,2))
        print(p4,vp=vplayout(3,1:2))
        print(p8,vp=vplayout(4,1:2))
        print(p7,vp=vplayout(5,1:2))
      }
    }
  }


}

