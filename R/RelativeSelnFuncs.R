## Weighted polynomial functions

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
  Poly=function(Catch,vars=c("lgth","nC","nT"),Quasi=F,wgt="AICc",All=TRUE) {
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

#' Fit model averaged 4th order polynomial to raw subsampled data
#' @description Fit model averaged 4th order polynomial to raw subsampled data
#'
#' @param Catch Matrix including columns "freqC", "SFC", "freqT" and "SFT"
#' @param fitTots Logical indicating whether to fit to totals (over hauls)
#'
#' @export
fitPoly=function(RawDf,fitTots=T,...) {
  RawDf=transform(RawDf,nC=freqC*SFC, nT=freqT*SFT)
  if(!fitTots) Df=RawDf[,c("lgth","nC","nT")]
     else Df=RawDf %>% group_by(lgth) %>% summarize(nC=sum(nC),nT=sum(nT))
  Df=round(as.data.frame(Df))
  Poly(Df,...)
}


#' Fit smooth function to catch share data
#' @description Fit spline to catch share
#'
#' @param Catch Matrix with data in first three columns
#' @param quasi Logical, whether to apply quasibinomial correction
#' @param bs Choice of smoother
#' @param k Dimension of the basis. k=3 in minimum for natural cubic spline.
#' @param m Order of the penalty
#' @param sp Supplied smoothing parameter
#'
#' @return List containing fitted model.
#' @export
FitGAM=function(Catch,quasi=F,bs="cr",k=5,m=NA,sp=NULL,null=F,rm.zeros=T) {
  lenname=colnames(Catch)[1]
  if(is.null(lenname)) lenname="lgth"
  formla=as.formula( paste0("y~s(",lenname,",bs=bs,k=k)") )
  #Catch$lgth=Catch[,1]
  Catch$n=Catch[,2]+Catch[,3]
  Catch$y=Catch[,3]/Catch$n
  Catch$y[is.na(Catch$y)]=0
  if(rm.zeros) Catch=Catch[Catch$n>0,]
  fam=ifelse(quasi,"quasibinomial","binomial")
  Gam.fit=gam(formla,family=fam,weights=n,sp=sp,data=Catch)
  #Gam.fit=gam(y~s(lgth,bs=bs,k=k),family=fam,weights=n,sp=sp,data=Catch)
  if(null) Gam.fit=gam(y~1,family=fam,weights=n,sp=sp,data=Catch)
  return(Gam.fit)
}


#' Fit logistic seln curve to paired-haul data
#' @description Fit logistic seln curve to paired-haul data to raw subsampled data
#'
#' @param Catch Matrix including columns "freqB", "SF.B", "freqG" and "SF.G"
#' @param fitTots Logical indicating whether to fit to totals (over hauls)
#'
#' @export
fitPH=function(RawDf,fitTots=T) {
  RawDf=transform(RawDf,nB=freqB*SF.B, nG=freqG*SF.G)
  if(!fitTots) Df=RawDf[,c("lgth","nB","nG")]
     else Df=RawDf %>% group_by(lgth) %>% summarize(nfine=sum(nB),nwide=sum(nG))
  Df=round(as.data.frame(Df))
  Fit=SELECT(Df,dtype="ph",print.out=F,x0=c(-5,0.5,0),penalty.func=NULL)
}

qAIC=function(fit,OD,correct=F) {
  Correction=0
  if(correct) {
    k=fit$rank
    Correction=2*k*(k+1)/(fit$df.resid-1)
  }
  qAIC=-2*logLik(fit)/OD+2*(attributes(logLik(fit))$df)+Correction
}

CalcOD=function(Df,yhat,suff.big=3) { #Assumes (propn) Df$y and Df$n exist
  y=Df$y; n=Df$n
  PearsonChisq=(y-yhat)^2/(yhat*(1-yhat)/n)
  Include=(n*yhat>suff.big & n*(1-yhat)>suff.big)
  nbig=sum(Include)
  PearsonStat=sum(PearsonChisq[Include])/(nbig-5)
  list(PearsonStat,nbig)
}

