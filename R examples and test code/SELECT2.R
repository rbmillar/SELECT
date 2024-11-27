
SELECT2=function(data,var.names,dtype="cc",stype="logistic",useTots=TRUE,
          q.names=NULL,Meshsize=NULL,x0=NULL,rel.power=NULL,penalty.func=NULL,
          verbose=T,control=list(maxit=10000,reltol=1e-8),Fit=T) {
  if(typeof(var.names)!="character")
    stop('SELECT errror message: \n Variable names must be character')
  if(!is.null(q.names) & typeof(q.names)!="character")
    stop('SELECT errror message: \n Scaling variable names must be character')
  dtype=substr(dtype,1,2)
  rtype=paste0(dtype,".",stype)
  r=propncurves(rtype) #Get propn catch curve function
  nGears=length(var.names)-1
  GearNames=paste0("n",1:nGears)
  Data=data[,var.names]
  colnames(Data)=c("lgth",GearNames)
  if(!is.null(q.names)) Data[,-1] = Data[,-1]/data[,q.names]
  if(useTots) Data=Data %>% group_by(lgth) %>%
           summarize(across(all_of(GearNames),sum)) %>% data.frame()
  Counts=Data[,-1]
  Const=mchoose(Counts,log=T) #Constant for log-likelihoods
  designTypes=c("relative","covered codend","paired haul","paired haul")
  design=designTypes[match(dtype,c("re","cc","ph","un"))]
  if(!(dtype %in% c("re","cc","ph","un")))
    stop('SELECT errror message: \ndesign must be one of "cc" (covered codend),
     "ph" (paired haul), or "re" (relative - primarily for gillnets and hooks)')
  rtype=paste0(dtype,".",stype)
  if(stype %in% c("logistic","richards") & is.null(Meshsize)) Meshsize=c(0,1)
  if(dtype=="un") Meshsize=c(1,1) #Arbitrary, since unrestricted
  if(sum(sort(Meshsize)==Meshsize)!=length(Meshsize))
    stop("Mesh size must be in ascending order")
  if(is.null(rel.power) | dtype %in% (c("cc","ph","un")))
    rel.power=rep(1,length(Meshsize))
  if(is.null(penalty.func)) penalty.func=function(theta){0.0}
  if(nGears!=length(Meshsize)) stop("Number of mesh sizes should be ",nGears)
  if(is.null(x0)) x0=StartValues(rtype,Data)
  SELECT.args=list(Data=Data,rtype=rtype,Meshsize=Meshsize,x0=x0,
                   rel.power=rel.power,penalty.func=penalty.func)
  #Calculate logliks at x0 and of saturated model
  ll.init=Const-nllhood(theta=x0,Data,Meshsize,r,rel.power,penalty.func)
  CountTotals=apply(Counts,1,sum,na.rm=TRUE)
  npos=sum(CountTotals>0)
  CountTotals=ifelse(CountTotals==0,Inf,CountTotals)
  CountPropns=Counts/CountTotals
  ll.fullfit=Const+sum(Counts[Counts>0]*log(CountPropns[Counts>0]),na.rm=TRUE)
  if(verbose) cat("Log-likelihood is",ll.init,"at x0=",round(x0,2),"\n")
  if(verbose) cat("Saturated log-likelihood is",ll.fullfit,"\n")
  if(!Fit) {
    if(verbose) cat("SELECT model fitted at x0 - no optimization: \n")
    control$maxit=0
    #Hessian is not evaluated at x0 when maxit=0, so set hessian=F
    fit=optim(x0,nllhood,Data=Data,Meshsize=Meshsize,r=r,rel.power=rel.power,
              penalty.func=penalty.func,hessian=F,control=control)
    fit$par=x0 #Quirk of optim with maxit=0 is that pars are set to zero
    fit$hessian=matrix(NA,length(x0),length(x0)) }
  if(Fit) {
    if(verbose)
      cat("Fitting SELECT model with",stype,"selection curves to",design,"data.\n")
    fit=optim(x0,nllhood,Data=Data,Meshsize=Meshsize,r=r,rel.power=rel.power,
            penalty.func=penalty.func,hessian=T,control=control) }
  fit$value=fit$value-Const
  fit=fit[names(fit) %in% c("counts", "message") == FALSE]
  Dev=2*(ll.fullfit+fit$value)
  if(verbose) {
    cat(paste0("Convergence code ",fit$convergence,
             ": Optimizer has ",ifelse(fit$convergence==0,"","*NOT*"),"converged\n"))
    cat("Pars=",fit$par,", Deviance=",Dev,", #len classes=",npos,"\n") }
  z=c(Call=match.call(),fit,deviance=Dev,SELECT.args,npos=npos,
         init.logLik=ll.init,full.logLik=ll.fullfit,logLik=-fit$value)
  class(z)="SELECT"
  return(invisible(z))
#  invisible(c(fit,deviance=Dev,rtype=rtype,rel.power=list(rel.power),
#              Meshsize=list(Meshsize),Data=list(Data)))
#  invisible(list(fit=fit,deviance=Dev,Call=SELECTcall,Data=Data,rtype=rtype,
#              rel.power=rel.power,Meshsize=Meshsize) )
}

#S2=SELECT2(c("lgth","n1","n2"),"ph",data=Catch)
#S1=SELECT(Catch[,1:3],"ph")

#S2R=SELECT2(c("lgth","n1","n2"),"ph",data=Catch,stype="richards")
#S1R=SELECT(Catch[,1:3],"ph",stype="richards")

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
#summary(Fit1)

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

