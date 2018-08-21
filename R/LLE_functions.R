#' @title predicttsir
#' @description function to predict incidence and susceptibles using the tsir model. This is different than
#' simulatetsir as you are inputting parameters as vectors. The output is a data frame I and S with mean and confidence intervals of predictions.
#' @param times The time vector to predict the model from. This assumes that the time step is equal to IP
#' @param births The birth vector (of length length(times) or a single element) where each element is the births in that given (52/IP) time step
#' @param beta The length(52/IP) beta vector of contact.
#' @param alpha A single numeric which acts as the homogeniety parameter.
#' @param S0 The starting initial condition for S. This should be greater than one, i.e. not a fraction.
#' @param I0 The starting initial condition for I. This should be greater than one, i.e. not a fraction.
#' @param nsim The number of simulations to perform.
#' @param stochastic A TRUE / FALSE argument where FALSE is the deterministic model, and TRUE is a negative binomial distribution.
#' @export
#' @examples
#' \dontrun{
#' require(kernlab)
#' require(ggplot2)
#' require(kernlab)
#' require(tsiR)
#' London <- twentymeas$London
#'
#' London <- subset(London, time > 1950)
#'
#' IP <- 2
#' ## first estimate paramters from the London data
#' parms <- estpars(data=London, IP=2, regtype='gaussian')
#'
## look at beta and alpha estimate
#' plotbeta(parms)
#'
#' ## now lets predict forward 20 years using the mean birth rate,
#' ## starting from rough initial conditions
#' births <- min(London$births)
#' times <- seq(1965,1985, by = 1/ (52/IP))
#' S0 <- parms$sbar
#' I0 <- 1e-5*mean(London$pop)
#'
#' pred <- predicttsir(times=times,births=births,
#'                     beta=parms$contact$beta,alpha=parms$alpha,
#'                     S0=S0,I0=I0,
#'                     nsim=50,stochastic=T)
#'
#' ## plot this prediction
#' ggplot(pred$I,aes(time,mean))+geom_line()+geom_ribbon(aes(ymin=low,ymax=high),alpha=0.3)
#'
#'
#'}

predicttsir <- function(times,births,beta,alpha,S0,I0,nsim,stochastic){

  I.mat <- S.mat <-  matrix(NA,length(times),nsim)

  alpha <- alpha[1]

  if(length(beta) < length(times)){
    beta <- rep(beta,length(times))[1:length(times)]
  }

  if(length(births) < length(times)){
    births <- rep(births,length(times))
  }

  for(n in 1:nsim){

    S <- I <- rep(NA,length(times))
    S[1] <- round(S0)
    I[1] <- round(I0)

    for(t in 2:length(times)){

      lambda <- min(S[t-1],unname(beta[t-1] * S[t-1] * (I[t-1])^alpha))

      if(stochastic){
        I[t] <- rnbinom(n=1,mu=lambda,size=I[t-1]+1e-10)
      }else{
        I[t] <- lambda
      }
      S[t] <- max(S[t-1] + births[t-1] - I[t],0)

    }
    I.mat[,n] <- I
    S.mat[,n] <- S
  }

  I.mat <- data.frame(I.mat)
  S.mat <- data.frame(S.mat)

  I.error <- apply(as.matrix(I.mat), 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

  I.mat$high <- I.error[2,]
  I.mat$low <- I.error[1,]

  S.error <- apply(as.matrix(S.mat), 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

  S.mat$high <- S.error[2,]
  S.mat$low <- S.error[1,]

  I.mat$mean <- rowMeans(I.mat[,1:nsim],na.rm=T)
  S.mat$mean <- rowMeans(S.mat[,1:nsim],na.rm=T)

  I.mat$time <- times
  S.mat$time <- times

  return(list('I'=I.mat,'S'=S.mat))
}




#' @title TSIR_LE
#' @description A function to calculate the Lyapunov Exponennt (LE) from the TSIR model
#' @param I The I output from the simulated or predicted TSIR model
#' @param S The S output from the simulated or predicted TSIR model
#' @param alpha The homogeneity parameter from the simulated or predicted TSIR model
#' @param beta The inferred contact rate, use beta = contact$beta where contact is an output from runtsir or simulatetsir
#' @param time The time vector from the data or simulated data
#' @param IP The generation interval of the pathogen (in weeks)
#'
#' @examples
#'
#' \dontrun{
#' require(kernlab)
#' require(ggplot2)
#' require(kernlab)
#' London <- twentymeas$London
#' ## just analyze the biennial portion of the data
#' London <- subset(London, time > 1950)
#'
#' ## define the interval to be 2 weeks
#' IP <- 2
#'
#'## first estimate paramters from the London data
#' parms <- estpars(data=London, IP=2, regtype='gaussian',family='poisson',link='log')
#'
#' ## look at beta and alpha estimate
#' plotbeta(parms)
#'
#' ## simulate the fitted parameters
#' sim <- simulatetsir(data=London,parms=parms,IP=2,method='deterministic',nsim=2)
#'
#' ## now lets predict forward 200 years using the mean birth rate,
#' ## starting from rough initial conditions
#' times <- seq(1965,2165, by = 1/ (52/IP))
#' births <- rep(mean(London$births),length(times))
#' S0 <- parms$sbar
#' I0 <- 1e-5*mean(London$pop)
#'
#' pred <- predicttsir(times=times,births=births,
#'                    beta=parms$contact$beta,alpha=parms$alpha,
#'                   S0=S0,I0=I0,
#'                   nsim=50,stochastic=T)
#'
#' ## take the last 10 years
#' pred <- lapply(pred, function(x)  tail(x, 52/IP * 20) )
#'
#' ## now compute the Lyapunov Exponent for the simulate and predicted model
#'
#' simLE <- TSIR_LE(
#' time=sim$res$time,
#' S=sim$simS$mean,
#' I=sim$res$mean,
#' alpha=sim$alpha,
#'   beta=sim$contact$beta,
#' IP=IP
#' )
#'
#' predLE <- TSIR_LE(
#' time=pred$I$time,
#' S=pred$S$X3,
#' I=pred$I$X3,
#' alpha=parms$alpha,
#' beta=parms$contact$beta,
#' IP=IP
#' )
#'
#' simLE$LE
#' predLE$LE
#'
#'
#' }

TSIR_LE = function(time, S, I, alpha, beta, IP) {
  IT <- length(I)
  j11 <- rep(NA, IT)
  j12 <- rep(NA, IT)
  j21 <- rep(NA, IT)
  j22 <- rep(NA, IT)
  # initial unit vector
  J <- matrix(c(1, 0), ncol = 1)
  # loop over the attractor
  for (i in 1:IT) {
    j11[i] <- 1 - beta[i %% (52/IP) + 1] * I[i]^alpha
    j12[i] <- - beta[i %% (52/IP) + 1] * S[i] * (I[i]^(alpha - 1) * alpha)
    j21[i] <- beta[i %% (52/IP) + 1] * I[i]^alpha
    j22[i] <- beta[i %% (52/IP) + 1] * S[i] * I[i]^(alpha - 1) * alpha
    J <- matrix(c(j11[i], j12[i], j21[i], j22[i]), ncol = 2, byrow = TRUE) %*% J
    if(is.nan(J[1])) stop(paste0('Matrix product has become NaN at index ',i))
  }
  res <- list(LE = log(norm(J))/IT, j11 = j11, j12 = j12, j21 = j21, j22 = j22,time=time,
              I = I, S = S, alpha = alpha, beta = beta)
  return(res)
}




#' @title TSIR_LLE
#' @description A function to calculate the Local Lyapunov Exponennt (LLE) from the TSIR model
#' @param LE The output of TSIR_LE to pass the Jacobian elements
#' @param m The window to sweep the time-varying Jacobian elements. Defaults to one.
#' @examples
#'
#'
#' \dontrun{
#' require(kernlab)
#' require(ggplot2)
#' require(kernlab)
#' London <- twentymeas$London
#' ## just analyze the biennial portion of the data
#' London <- subset(London, time > 1950)
#'
#' ## define the interval to be 2 weeks
#' IP <- 2
#'
#'## first estimate paramters from the London data
#' parms <- estpars(data=London, IP=2, regtype='gaussian',family='poisson',link='log')
#'
#' ## look at beta and alpha estimate
#' plotbeta(parms)
#'
#' ## simulate the fitted parameters
#' sim <- simulatetsir(data=London,parms=parms,IP=2,method='deterministic',nsim=2)
#'
#' ## now lets predict forward 200 years using the mean birth rate,
#' ## starting from rough initial conditions
#' times <- seq(1965,2165, by = 1/ (52/IP))
#' births <- rep(mean(London$births),length(times))
#' S0 <- parms$sbar
#' I0 <- 1e-5*mean(London$pop)
#'
#' pred <- predicttsir(times=times,births=births,
#'                    beta=parms$contact$beta,alpha=parms$alpha,
#'                   S0=S0,I0=I0,
#'                   nsim=50,stochastic=T)
#'
#' ## take the last 10 years
#' pred <- lapply(pred, function(x)  tail(x, 52/IP * 20) )
#'
#' ## now compute the Lyapunov Exponent for the simulate and predicted model
#'
#' simLE <- TSIR_LE(
#' time=sim$res$time,
#' S=sim$simS$mean,
#' I=sim$res$mean,
#' alpha=sim$alpha,
#'   beta=sim$contact$beta,
#' IP=IP
#' )
#'
#' predLE <- TSIR_LE(
#' time=pred$I$time,
#' S=pred$S$X3,
#' I=pred$I$X3,
#' alpha=parms$alpha,
#' beta=parms$contact$beta,
#' IP=IP
#' )
#'
#' simLE$LE
#' predLE$LE
#'
#' simLLE <- TSIR_LLE(simLE)
#' predLLE <- TSIR_LLE(predLE)
#'
#' plotLLE(simLLE)
#' plotLLE(predLLE)
#' }
#'


TSIR_LLE = function(LE, m = 1) {
  LLE <- rep(NA, length(LE$I))
  for (i in 1:(length(LE$I) - m)) {
    J <- matrix(c(1, 0, 0, 1), ncol = 2)
    for (k in 0:(m - 1)) {
      J <- matrix(c(LE$j11[(i + k)], LE$j12[(i + k)], LE$j21[(i + k)], LE$j22[(i + k)]), ncol = 2, byrow = TRUE) %*% J
    }
    LLE[i] = log(max(abs(eigen(J)$values)))/m
  }
  res <- list(LLE = LLE, time=LE$time, I = LE$I, S = LE$S)
  return(res)
}


#' @title plotLLE
#' @description Function to plot the Local Lyapunov Exponents. The output is of class ggplot2 so you can add standard
#' ggplot2 options to it if desired.
#' @param LLE The output from TSIR_LLE
#' @export
#' @examples
#'
#'
#' \dontrun{
#' require(kernlab)
#' require(ggplot2)
#' require(kernlab)
#' London <- twentymeas$London
#' ## just analyze the biennial portion of the data
#' London <- subset(London, time > 1950)
#'
#' ## define the interval to be 2 weeks
#' IP <- 2
#'
#'## first estimate paramters from the London data
#' parms <- estpars(data=London, IP=2, regtype='gaussian',family='poisson',link='log')
#'
#' ## look at beta and alpha estimate
#' plotbeta(parms)
#'
#' ## simulate the fitted parameters
#' sim <- simulatetsir(data=London,parms=parms,IP=2,method='deterministic',nsim=2)
#'
#' ## now lets predict forward 200 years using the mean birth rate,
#' ## starting from rough initial conditions
#' times <- seq(1965,2165, by = 1/ (52/IP))
#' births <- rep(mean(London$births),length(times))
#' S0 <- parms$sbar
#' I0 <- 1e-5*mean(London$pop)
#'
#' pred <- predicttsir(times=times,births=births,
#'                    beta=parms$contact$beta,alpha=parms$alpha,
#'                   S0=S0,I0=I0,
#'                   nsim=50,stochastic=T)
#'
#' ## take the last 10 years
#' pred <- lapply(pred, function(x)  tail(x, 52/IP * 20) )
#'
#' ## now compute the Lyapunov Exponent for the simulate and predicted model
#'
#' simLE <- TSIR_LE(
#' time=sim$res$time,
#' S=sim$simS$mean,
#' I=sim$res$mean,
#' alpha=sim$alpha,
#'   beta=sim$contact$beta,
#' IP=IP
#' )
#'
#' predLE <- TSIR_LE(
#' time=pred$I$time,
#' S=pred$S$X3,
#' I=pred$I$X3,
#' alpha=parms$alpha,
#' beta=parms$contact$beta,
#' IP=IP
#' )
#'
#' simLE$LE
#' predLE$LE
#'
#' simLLE <- TSIR_LLE(simLE)
#' predLLE <- TSIR_LLE(predLE)
#'
#' plotLLE(simLLE)
#' plotLLE(predLLE)
#' }
#'

plotLLE <- function(LLE){

  local.ly <- LLE$LLE

  pm <- LLE$LLE > 0

  S <- LLE$S

  I <- LLE$I

  time <- LLE$time

  LLdf <- data.frame(time,S,I,'LLEsign'=pm,'LLE'=local.ly)

  LLdf <- head(LLdf,-1)

  ggplot(LLdf,aes(time,I,colour=LLE,group=1))+geom_line()

  p1 <- ggplot(LLdf,aes_string('S','I',size='abs(LLE)',shape='LLEsign',colour='LLEsign'))+geom_point(alpha=0.7)+ theme_bw()+
    theme(text=element_text(size=12))+
    scale_size_continuous(limits = c(0,1),range = c(2,10),name = 'abs(LLE)')+
    scale_shape_manual(values=c(1,19),labels = c("LLE < 0", "LLE > 0"),name='sign(LLE)')+
    # labs(size="abs(LLE)",
    #      col='sign(LLE)', shape='sign(LLE)')+
    #scale_color_manual(values=c('grey2','grey'))+
    scale_color_manual(values=c('orangered3','dodgerblue'),labels = c("LLE < 0", "LLE > 0"),name='sign(LLE)')
  #guides(shape=FALSE,size=FALSE,colour=FALSE)

  return(p1)
}


