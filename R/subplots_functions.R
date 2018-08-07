#' @title plotregression
#'
#' @description Plots the cumulative cases - cumulative births data and regression fit
#' @param dat the list produced from the runtsir, mcmctsir, and simulatetsir function.
plotregression <- function(dat){

  regdf <- NULL
  regdf$X <- dat$X
  regdf$Y <- dat$Y
  regdf$Yhat <- dat$Yhat
  regdf$time <- dat$res$time
  regdf <- as.data.frame(regdf)
  meltdf <- melt(regdf, id.vars = "time")
  p1 <- ggplot(meltdf, aes_string("time", "value", col = "variable")) +
    geom_line(size = 2) + theme_bw()

  print(p1)
}

#' @title plotrho
#'
#' @description Plots the inferred reporting rate, rho
#' @param dat the list produced from the runtsir, mcmctsir, and simulatetsir function.
plotrho <- function(dat){

  rhodf <- NULL
  rhodf$time <- dat$res$time
  rhodf$rho <- 1/dat$rho
  rhodf <- as.data.frame(rhodf)
  p2 <- ggplot(rhodf, aes_string("time", "rho")) + geom_line(size = 2) +
    theme_bw() + ggtitle(bquote(bar(rho) == .(signif(mean(1/dat$rho),
                                                     2)))) + ylab(bquote(1/rho))

  print(p2)
}

#' @title plotsbar
#'
#' @description Plots the profile log likelihood calculation for inferred sbar
#' @param dat the list produced from the runtsir, mcmctsir, and simulatetsir function.
plotsbar <- function(dat){

  loglikdf <- NULL
  loglikdf$sbar <- dat$Smean
  loglikdf$loglik <- dat$loglik
  loglikdf <- as.data.frame(loglikdf)
  p9 <- ggplot(loglikdf, aes_string("sbar", "loglik")) +
    geom_line() + geom_point() + theme_bw() + geom_vline(xintercept = dat$sbar,  linetype = "longdash") +
    ggtitle(bquote(bar(S) ==    .(signif(dat$sbar, 2)) ~ "," ~ .(signif(dat$sbar/mean(dat$pop) * 100, 2)) ~ "%")) + xlab(bquote(bar(S)))

  print(p9)
}

#' @title plotbeta
#'
#' @description Plots the inferred beta with confidence intervals (when they can be calculated)
#' @param dat the list produced from the runtsir, mcmctsir, and simulatetsir function.
plotbeta <- function(dat){
  # betadf <- NULL
  # betadf$time <- seq(1, length(dat$beta), 1)
  # betadf$beta <- dat$beta

  betadf <- dat$contact

  betadf <- as.data.frame(betadf)
  p4 <- ggplot(betadf, aes_string("time", "beta")) + geom_line(size = 2) +
    theme_bw() + ggtitle(bquote(bar(beta) == .(signif(mean(dat$beta),
                                                      2)) ~ "," ~ alpha == .(signif(dat$alpha, 3)))) +
    ylab(bquote(beta))
  if ("contact" %in% names(dat)) {
    p4 <- ggplot(dat$contact, aes_string("time", "beta")) +
      geom_line(size = 2) + geom_ribbon(ymin = dat$contact$betalow,
                                        ymax = dat$contact$betahigh, alpha = 0.5, col = "dodgerblue",
                                        fill = "dodgerblue") + ylim(c(min(dat$contact$betalow),
                                                                      max(dat$contact$betahigh))) + theme_bw() + ggtitle(bquote(bar(beta) ==
                                                                                                                                  .(signif(mean(dat$beta), 2)) ~ "," ~ alpha ==
                                                                                                                                  .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
    if (sum(sum(is.na(dat$contact))) > 0) {
      p4 <- ggplot(betadf, aes_string("time", "beta")) +
        geom_line(size = 2) + theme_bw() + ggtitle(bquote(bar(beta) ==
                                                            .(signif(mean(dat$beta), 2)) ~ "," ~ alpha ==
                                                            .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
    }
  }
  p4 <- p4 + xlab(sprintf("time mod %g", nrow(dat$contact)))


  print(p4)
}

#' @title plotforward
#'
#' @description Plots the forward simulation from the TSIR model
#' @param dat the list produced from the runtsir, mcmctsir, and simulatetsir function.
#' @param inverse a TRUE or FALSE option to plot the forward simulate negative (TRUE) or positive (FALSE). Defaults to FALSE.
plotforward <- function(dat,inverse = F){
  drops <- c("mean", "sd", "error", "cases", "time")
  sim.only <- dat$res[, !(names(dat$res) %in% drops)]
  n <- ncol(sim.only)
  error <- qt(0.975, df = n - 1) * dat$res$sd/sqrt(n)
  dat$res$error <- error
  eb <- aes(ymax = mean + error, ymin = mean - error)
  p6 <- ggplot(data = dat$res, aes_string("time")) + theme(legend.position = "none") +
    geom_line(aes_string(y = "cases"), colour = "dodgerblue",
              size = 1) + xlab("year") + ylab("cases") + geom_line(aes_string(y = "mean"),
                                                                   colour = "orangered4", size = 1) + geom_ribbon(eb,
                                                                                                                  alpha = 0.3) + theme_bw()
  inversecases <- dat$res
  inversecases$cases <- -dat$res$cases
  p7 <- ggplot(data = inversecases, aes_string("time")) +
    theme(legend.position = "none") + geom_line(aes_string(y = "cases"),
                                                colour = "dodgerblue", size = 1) + xlab("time") +
    ylab("cases") + geom_line(aes_string(y = "mean"),
                              colour = "orangered4", size = 1) + geom_ribbon(eb,
                                                                             alpha = 0.3) + theme_bw()
  if(inverse){
    print(p7)
  }else{
    print(p6)
  }
}
