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
  Wk=ModelCheck(obj,print.out=F,plots=c(F,F))
  FitStats=Wk$stats
  GOFStats=Wk$gof
  llhoods=FitStats[,c("null.l","model.l","full.l","AIC")]
  cat("\nLog-likelihoods and AIC:\n")
  print(llhoods)
  cat("\nGoodness of fit statistics\n")
  GOF=GOFStats[,c("Deviance","Pearson.chisq","dof","Deviance.CF","Pearson.CF")]
  print(GOF)
  cat("\nEstimates\n")
  print(Estimates(obj))
  invisible(list(Call=obj$Call,llhoods=llhoods,GOF=GOF))
}
#summary(Fit1)

#' @export
plot.SELECT=function(obj,plotlens=NULL,npts=101,...) {
  if(is.null(plotlens)) {
    lgth=obj$Data[,1]
    plotlens=seq(min(lgth),max(lgth),length=npts) }
  PlotCurves(obj,plotlens,...)
}

#' @export
predict.SELECT=function(obj,lens=NULL,...) {
  if(!is.null(lens) & !is.vector(lens))
    stop("Error: lens must be a vector of lengths")
  PlotCurves(obj,plotlens=lens,plot.out=F,...) }

#' @export
fitted.SELECT=function(obj,...) {
  Fit=PlotCurves(obj,plot.out=F,...) 
  Fit }

#' @export
logLik.SELECT=function(obj,type="SELECT") {
  logLik=obj$logLik
  attr(logLik,"nobs")=obj$npos
  class(logLik)="logLik"
  if(type=="SELECT") attr(logLik,"df")=length(obj$par)
  else
  { O=as.matrix(obj$Data[,-1])
    E=ModelCheck(obj,print.out=F,plots=c(F,F))$fit
    attr(logLik,"df")=length(obj$par)+nrow(O)
    logLik=sum(dpois(as.vector(O),as.vector(E),log=T))
    names(logLik)="Poisson.logLik" }
  logLik }

#' @export
deviance.SELECT=function(obj) obj$deviance

#' Calculate AIC
#' @description Calculates AIC of SELECT object. By default it provides the
#' binomial SELECT AIC. For comparison with models fitted by glm of gam use
#' type="poisson" to get the Poisson AIC.
#'
#' @param obj Fitted SELECT model
#' @param type Use type="poisson" to get Poisson AIC. This is calculated from the
#' matrices of observed and fitted counts and does not condition on row totals.
#' @export
AIC.SELECT=function(obj,type="SELECT") {
  if(type=="SELECT") AIC=ModelCheck(obj,print.out=F,plots=c(F,F))$stats[1,"AIC"]
  else
    { O=as.matrix(obj$Data[,-1])
      npar=length(obj$par)+nrow(O)
      E=ModelCheck(obj,print.out=F,plots=c(F,F))$fit
      AIC=-2*sum(dpois(as.vector(O),as.vector(E),log=T)) + 2*npar
      names(AIC)="Poisson.AIC" }
  AIC }

#' Calculate BIC
#' @description Calculates BIC of SELECT object. By default it provides the
#' binomial SELECT BIC. For comparison with models fitted by glm of gam use
#' type="poisson" to get the Poisson BIC.
#'
#' @param obj Fitted SELECT model
#' @param type Use type="poisson" to get Poisson BIC. This is calculated from the
#' matrices of observed and fitted counts and does not condition on row totals.
#' @export
BIC.SELECT=function(obj,type="SELECT") {
  if(type=="SELECT") BIC=ModelCheck(obj,print.out=F,plots=c(F,F))$stats[1,"BIC"]
  else
  { O=as.matrix(obj$Data[,-1])
  nobs=nrow(O)*ncol(O) #Note: All-zero rows are data, so not excluded
  npar=length(obj$par)+nrow(O)
  E=ModelCheck(obj,print.out=F,plots=c(F,F))$fit
  BIC=-2*sum(dpois(as.vector(O),as.vector(E),log=T)) + log(nobs)*npar
  names(BIC)="Poisson.BIC" }
  BIC }


