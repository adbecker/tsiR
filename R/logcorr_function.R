#' @title logcorr
#'
#' @description Plot the correlation of the true data against the fitted resimulated data.
#'
#' @param sim The dataframe or list produced by the 'runtsir' function.

logcorr <- function(sim){
  
  if(class(sim) == "list"){
    sim <- sim$res
  }
  if(class(sim) == "data.frame"){
    sim <- sim
  }
  
  obs <- log(sim$cases+1)
  pred <- log(sim$mean+1)
  fit <- lm(pred ~ obs)
  rsquared <- signif(summary(fit)$adj.r.squared, 2)
  pval <- signif(summary(fit)$coef[2,4], 2)
  c1 <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() + stat_smooth(method = "lm", col = "black") +theme_bw()+
    ggtitle(bquote(Adjusted~' '~R^{2}==.(rsquared)))+xlab('log(observed)')+ylab('log(predicted)')
  return(c1)
}

