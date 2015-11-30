ggreg <- function(fit){
  rsquared <- signif(summary(fit)$adj.r.squared, 2)
  pval <- signif(summary(fit)$coef[2,4], 2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() + stat_smooth(method = "lm", col = "black") +theme_bw()+
    ggtitle(bquote(Adjusted~' '~R^{2}==.(rsquared)))
}