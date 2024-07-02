## Functions to fit SELECT models, including spline and model-averaged poly
# SELECT, splineSELECT & polySELECT

#===============================================================================
#' Fit selection curves to fishing gears using the SELECT method.
#' @description SELECT is an acronym for the Share Each LEngths Catch Total method
#' of Millar (1992). SELECT is a general statistical framework to  estimate the
#' selectivity of trawls, hooks and gillnets from experimental data.
#' @references Millar, R. B. (1992).
#' Estimating the size-selectivity of fishing gear by conditioning on the
#' total catch. Journal of the American Statistical Association. 87: 962-968.
#' @export
SELECT=function(data,var.names,dtype="cc",stype="logistic",useTots=TRUE,
                 q.names=NULL,Meshsize=NULL,x0=NULL,rel.power=NULL,penalty.func=NULL,
                 verbose=NULL,control=list(maxit=10000,reltol=1e-8),Fit=TRUE) {
  #SELECT.args=as.list(environment()) %>% discard(is.null)
  if(typeof(var.names)!="character")
    stop('SELECT errror message: \n Variable names must be character')
  if(!is.null(q.names) & typeof(q.names)!="character")
    stop('SELECT errror message: \n Scaling variable names must be character')
  dtype=substr(dtype,1,2)
  rtype=paste0(dtype,".",stype)
  r=propncurves(rtype) #Get propn catch curve function
  nGears=length(var.names)-1
  #Data=data[,var.names]
  #GearNames=paste0("n",1:nGears)
  #colnames(Data)=c("lgth",GearNames)
  #if(!is.null(q.names)) Data[,-1] = Data[,-1]/data[,q.names]
  #if(useTots) Data=Data %>% group_by(lgth) %>%
  #  summarize(across(all_of(GearNames),sum)) %>% data.frame()
  Data=Raw2Tots(data,var.names,q.names=q.names,useTots=useTots)
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
  SELECT.args=list(data=data,var.names=var.names,Data=Data,rtype=rtype,
        Meshsize=Meshsize,x0=x0,rel.power=rel.power,penalty.func=penalty.func)
  #Calculate logliks at x0 and of saturated model
  ll.init=Const-nllhood(theta=x0,Data,Meshsize,r,rel.power,penalty.func)
  CountTotals=apply(Counts,1,sum,na.rm=TRUE)
  npos=sum(CountTotals>0)
  CountTotals=ifelse(CountTotals==0,Inf,CountTotals)
  CountPropns=Counts/CountTotals
  ll.fullfit=Const+sum(Counts[Counts>0]*log(CountPropns[Counts>0]),na.rm=TRUE)
  if(is.null(verbose)) cat("Log-likelihood is",ll.init,"at x0=",round(x0,2),"\n")
  if(is.null(verbose)) cat("Saturated log-likelihood is",ll.fullfit,"\n")
  if(!Fit) {
    if(is.null(verbose)) cat("SELECT model fitted at x0 - no optimization: \n")
    control$maxit=0
    #Hessian is not evaluated at x0 when maxit=0, so set hessian=F
    fit=optim(x0,nllhood,Data=Data,Meshsize=Meshsize,r=r,rel.power=rel.power,
              penalty.func=penalty.func,hessian=F,control=control)
    fit$par=x0 #Quirk of optim with maxit=0 is that pars are set to zero
    fit$hessian=matrix(NA,length(x0),length(x0)) }
  if(Fit) {
    if(is.null(verbose))
      cat("Fitting SELECT model with",stype,"selection curves to",design,"data.\n")
    fit=optim(x0,nllhood,Data=Data,Meshsize=Meshsize,r=r,rel.power=rel.power,
              penalty.func=penalty.func,hessian=T,control=control) }
  fit$value=fit$value-Const
  fit=fit[names(fit) %in% c("counts", "message") == FALSE]
  Dev=2*(ll.fullfit+fit$value)
  if(is.null(verbose)) {
    cat(paste0("Convergence code ",fit$convergence,
               ": Optimizer has ",ifelse(fit$convergence==0,"","*NOT*"),"converged\n"))
    cat("Pars=",fit$par,", Deviance=",Dev,", #len classes=",npos,"\n") }
  z=c(Call=match.call(),fit,SELECT.args,npos=npos,deviance=Dev,
      init.logLik=ll.init,full.logLik=ll.fullfit,logLik=-fit$value)
  class(z)="SELECT"
  return(invisible(z))
}


#===============================================================================
#' Fit smooth function to catch share data
#' @description Fit spline to catch share - description needs updating
#'
#' @param Catch Matrix with data in first three columns
#' @param quasi Logical, whether to apply quasibinomial correction
#' @param bs Choice of smoother
#' @param k Dimension of the basis. k=3 in minimum for natural cubic spline.
#' @param m Order of the penalty
#' @param sp Supplied smoothing parameter
#'
#' @return List containing fitted model.
#'
#'
#' @export
SplineSELECT=function(data,var.names,q.names=NULL,useTots=TRUE,
                 bs="cr",k=5,sp=NULL,quasi=FALSE,rm.zeros=TRUE) {
  if(typeof(var.names)!="character")
    stop('SELECT errror: \n Variable names must be character')
  if(!is.null(q.names) & typeof(q.names)!="character")
    stop('SELECT errror \n Sampling fraction variable names must be character')
  Data=Raw2Tots(data,var.names,q.names=q.names,useTots=useTots)
  vn=var.names
  formla=as.formula( paste0("cbind(",vn[3],",",vn[2],")","~s(",vn[1],",bs=bs,k=k)") )
  if(rm.zeros) Data=Data[Data[,2]+Data[,3]>0,]
  fam=ifelse(quasi,"quasibinomial","binomial")
  GamFit=gam(formla,family=fam,sp=sp,data=Data)
  return(GamFit)
}

#===============================================================================
## Weighted polynomial functions - needs to be updated (20/5/2024)

#' Fit model averaged 4th order polynomial
#' @description Fit model averaged 4th order polynomial
#'
#' @param Catch Matrix with data in first three columns
#' @param vars Column variable names if not vars=c("lgth","nC","nT")
#' @param plens If specified, lengths at which to calculated fitted values.
#' @param Quasi Logical, whether to apply quasibinomial correction
#' @param wgt Weighting criterion. If not "AICc" then AIC is used
#' @param All If not TRUE then only full models from null to 4th order are used.
#'
#' @return List containing fitted values, overdispersion, wgt and dredge averaged fit.
#' @export
PolySELECT=function(Catch,vars=c("lgth","nC","nT"),Quasi=F,wgt="AICc",All=TRUE) {
  if(wgt!="AICc") wgt="AIC" #Only two possibilities for now
  if(All) Sub=expression(eval(TRUE))
  if(!All) Sub=expression( dc(I(x^0),x,I(x^2),I(x^3),I(x^4)) )
  Catch$x=Catch[,vars[1]]
  Catch$n=Catch[,vars[2]]+Catch[,vars[3]]
  Catch$y=Catch[,vars[3]]/Catch$n
  Catch=Catch[!is.na(Catch$y),]
  od=1
  fit4=glm(y~I(x^0)+x+I(x^2)+I(x^3)+I(x^4)-1,family=binomial,weights=n,data=Catch)
  if(!Quasi) {  fits=dredge(fit4 ,rank=wgt, subset=Sub) }
  if(Quasi) {
    qfit4=glm(y~I(x^0)+x+I(x^2)+I(x^3)+I(x^4)-1,family=quasibinomial,weights=n,
              data=Catch)
    od=summary(qfit4)$dispersion
    #od=CalcOD(Catch,fitted(fit4),1)[[1]]
    if(wgt=="AIC")
      fits=dredge(fit4, rank=function(fit) qAIC(fit,OD=od),subset=Sub)
    else fits=dredge(fit4, rank=function(fit) qAIC(fit,OD=od,correct=T),subset=Sub)
  }
  avg.fit=model.avg(fits,fit = T)
  list(avg.fit=avg.fit,od=od,wgt,fits)
}


#===============================================================================
## Log-likelihoods from full, null and fixed models

#' Evaluate log-likelihoods from full, null and fixed models
#' @export
LOGLIKS=function(data,var.names,q.names=NULL,useTots=TRUE,fixed=NULL) {
  if(typeof(var.names)!="character")
    stop('SELECT errror message: \n Variable names must be character')
  if(!is.null(q.names) & typeof(q.names)!="character")
    stop('SELECT errror message: \n Scaling variable names must be character')
  Data=Raw2Tots(data,var.names,q.names=q.names,useTots=useTots)
  Counts=Data[,-1]
  Const=mchoose(Counts,log=T) #Constant for log-likelihoods
  CountTotals=apply(Counts,1,sum,na.rm=TRUE)
  CountTotals=ifelse(CountTotals==0,Inf,CountTotals)
  CountPropns=Counts/CountTotals
  #Full fit llhood
  ll.fullfit=Const+sum(Counts[Counts>0]*log(CountPropns[Counts>0]),na.rm=TRUE)
  #Null fit llhood
  GearTotals=apply(Counts,2,sum,na.rm=TRUE)
  GearPropns=GearTotals/sum(GearTotals)
  ll.nullfit=Const+sum(GearTotals*log(GearPropns),na.rm=TRUE)
  #Null fit with fixed power llhood
  nGears=length(var.names)-1
  if(is.null(fixed)) fixed=rep(1,nGears)/nGears
  if(length(fixed)!=nGears | sum(fixed)!=1) {
    cat("\nError: fixed must be length",nGears,"and sum to unity\n")
    ll.nullfix=NA }
  else
    ll.nullfix=Const+sum(GearTotals*log(fixed),na.rm=TRUE)
  list(full=ll.fullfit,null=ll.nullfit,fixed=ll.nullfix)
}

