##' @title print.summary.mgee2
##' @param x the summary results
##' @param ... Other parameters
##' @return a table of summary statistics 
##' @export

print.summary.mgee2 = function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat('Summary table of the estimation\n')
  printCoefmat(x$TAB, digits=4)
}
