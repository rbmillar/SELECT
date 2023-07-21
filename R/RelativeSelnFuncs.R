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
  preds=predict(avg.fit,data.frame(x=lenseq), type="response")
  list(preds=preds,od=od,wgt,fits)
}

fitPoly=function(RawDf,fitTots=T,...) {
  RawDf=transform(RawDf,nC=freqC*SFC, nT=freqT*SFT)
  if(!fitTots) Df=RawDf[,c("lgth","nC","nT")]
     else Df=RawDf %>% group_by(lgth) %>% summarize(nC=sum(nC),nT=sum(nT))
  Df=round(as.data.frame(Df))
  Fit=Poly(Df,...)
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