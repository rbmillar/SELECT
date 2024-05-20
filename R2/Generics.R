print.SELECT=function(obj) {
  cat("Estimated parameters:",obj$par)
  cat("\nlog-likelihood:",obj$logLik)
  cat("\nDeviance:",obj$deviance) }

coef.SELECT=function(obj) Estimates(obj)[,"par"]


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

plot.SELECT=function(obj) ModelCheck(obj,print.out=F)

predict.SELECT=function(obj,newdata=NULL) {
  if(!is.null(newdata) & !is.vector(newdata))
     stop("Error: newdata must be a vector of lengths")
  PlotCurves(obj,plotlens=newdata,plot.out=F) }

logLik.SELECT=function(obj) {
  logLik=obj$logLik
  attr(logLik,"df")=length(obj$par)
  attr(logLik,"nobs")=obj$npos
  class(logLik)="logLik"
  logLik }

deviance.SELECT=function(obj) obj$deviance

