#' @title maxthreshold
#' @description A function used to optimize the threshold parameter to give the best fit to the data. Optimizes the fit based on R squared.
#' @param data The time, cases, births, pop data frame.
#' @param nsim The number of simulations to do.
#' @param IP The infectious period, which should the time step of the data.
#' @param method The forward simulation method used, i.e. deterministic, negbin, pois.
#' @param parms The estimated parameters from estpars or mcmcestpars.
#' @param thresholdmin The minimum number of cases to be considered an outbreak.
#' @param thresholdmax The max number of cases to be considered an outbreak.
#' @param printon A T/F statement to print the progress.
#' @param inits.fit Whether or not to fit initial conditions as well. Defaults to FALSE here. This parameter is more necessary in more chaotic locations.
#'
#' @examples
#' require(kernlab)
#' Mold <- twentymeas[["Mold"]]
#' plotdata(Mold)
#' \dontrun{
#'parms <- estpars(data=Mold,alpha=0.97)
#'tau <- maxthreshold(data=Mold,parms=parms,
#'thresholdmin=8,thresholdmax=12,inits.fit=FALSE)
#'res <- simulatetsir(data=Mold,parms=parms,
#'epidemics='break',threshold=tau,method='negbin',inits.fit=FALSE)
#' plotres(res)
#'}
maxthreshold <- function(data,nsim=2,IP=2,method='deterministic',inits.fit=FALSE,
                         parms,thresholdmin=2,thresholdmax=20,printon=FALSE){

  threshvec <- seq(thresholdmin,thresholdmax,1)
  rvec <- rep(0,length(threshvec))
  for(it in 1:length(rvec)){

    if(printon == T){
      print(sprintf('trying threshold=%d',threshvec[it]))
    }
    res <- simulatetsir(data, nsim = nsim, IP=IP,
                        parms=parms,
                        method=method,
                        inits.fit=inits.fit,
                        epidemics='break', pred ='forward',
                        threshold=threshvec[it],
                        add.noise.sd = 0, mul.noise.sd = 0)

    obs <- res$res$cases
    pred <- res$res$mean

    fit <- lm(pred ~ obs)
    r.squared <- summary(fit)$r.squared
    rvec[it] <- r.squared

    if(printon == T){
      print(sprintf('resulting Rsquared is %f',r.squared))
    }
  }

  threshold <- threshvec[which.max(rvec)]
  return(threshold)


}
