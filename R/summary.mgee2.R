#' summary.mgee2
#' @export
#' @param object The fitted model
#' @param ... Other parameters
#' summary function for mgee2 method output
summary.mgee2 <- function(object, ...){
  if(object$convergence==FALSE){
    cat("\nError: Model does not converge.\n")
    cat("\nCall:\n",
        paste(deparse(object$call), sep="\n", collapse = "\n"), "\n\n", sep="")
  }else{
    Est <- c(unname(object$beta),unname(object$alpha))
    Std.Err <- sqrt(diag(unname(object$variance)))
    z <- Est/Std.Err
    p <- 2*(1-pnorm(abs(z)))
    
    TAB <- list()
    TAB <- cbind(Est, Std.Err, z, p)
    colnames(TAB) = c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    rownames(TAB) = c(names(object$beta), names(object$alpha))

    results = list(call=object$call,TAB=TAB)
    class(results) = "summary.mgee2"
    return(results)
  }
}