##' @export
##' @import ggplot2
##' @importFrom graphics plot
##' @title plot_model
##' @description This function gives plot of the odds ratio or shows the iteration for convergence.
##' @param x results from the fitted model.
##' @param conv defulated for odds ratio plot, otherwise show the iteration plot. 
##' @return plot odds ratio with CIs or plot of the iterations.
##' @examples 
##'  beta=c(0.1,0.2,0.3)
##'  alpha=c(0.4,0.5)
##'  variance=c(0.8,0.5,0.7,0.3,0.4)
##'  x=list(beta,alpha,variance)
##'  names(x)=c("beta","alpha","variance")
##'  plot_model(x)

plot_model <- function(x, conv = FALSE){
  beta_length = length(x$beta)
  alpha_length = length(x$alpha)
  
  var = diag(x$variance)
  beta_var = var[1:beta_length]
  alpha_var = var[(beta_length+1):length(var)]
  
  odd_beta = exp(x$beta)
  odd_alpha = exp(x$alpha)
  
  zscore = qnorm(1-.05/2)
  beta_ci_lower = exp(x$beta-zscore*sqrt(beta_var))
  beta_ci_upper = exp(x$beta+zscore*sqrt(beta_var))
  alpha_ci_lower = exp(x$alpha-zscore*sqrt(alpha_var))
  alpha_ci_upper = exp(x$alpha+zscore*sqrt(alpha_var))
  
  df <- data.frame(
    x=seq(1:(beta_length+alpha_length)),
    y=c(odd_beta,odd_alpha),
    lower=c(beta_ci_lower,alpha_ci_lower),
    upper=c(beta_ci_upper,alpha_ci_upper))
  p <- ggplot(df, aes(x = x, y=y)) + labs(x = "order of estimators", y = "exp(point estimates)") + geom_point(colour="red")
  p= p + theme_bw()
  p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p = p + geom_linerange(aes(ymin = lower, ymax = upper),colour="blue")
  p = p + coord_flip()
  print(p)
  
  if (conv==TRUE){
    len <- length(x$differ)
    plot(0:(len-1), x$differ, xlab="Number of iterations", ylab="Difference")
  }
}