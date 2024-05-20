## Generic utility and extractor functions
# print, coef, summary, plot, predict, logLik & deviance

#' @export
print.SELECT=function(obj) {
  cat("Estimated parameters:",obj$par)
  cat("\nlog-likelihood:",obj$logLik)
  cat("\nDeviance:",obj$deviance) }

#' @export
coef.SELECT=function(obj) Estimates(obj)[,"par"]

#' @export
summary.SELECT=function(obj) {
  #Add AIC&BIC and table for par ests and sds
  cat("Call:",deparse(obj$Call),"\n")
  Estimates(obj)
  summStats=ModelCheck(obj,print.out=F,plots=c(F,F))$stats
  llhoods=summStats[,c("null.l","model.l","full.l")]
  cat("\nLog-likelihoods:\n")
  print(llhoods)
  cat("\nFit statistics\n")
  GOF=summStats[,c("Deviance","Pearson.chisq","dof","Deviance.CF","Pearson.CF")]
  print(GOF)
  cat("\nEstimates\n")
  print(Estimates(obj))
  invisible(list(Call=obj$Call,llhoods=llhoods,GOF=GOF))
}
#summary(Fit1)

#' @export
plot.SELECT=function(obj) ModelCheck(obj,print.out=F)

#' @export
predict.SELECT=function(obj,newdata=NULL) {
  if(!is.null(newdata) & !is.vector(newdata))
    stop("Error: newdata must be a vector of lengths")
  PlotCurves(obj,plotlens=newdata,plot.out=F) }

#' @export
logLik.SELECT=function(obj) {
  logLik=obj$logLik
  attr(logLik,"df")=length(obj$par)
  attr(logLik,"nobs")=obj$npos
  class(logLik)="logLik"
  logLik }

#' @export
deviance.SELECT=function(obj) obj$deviance

