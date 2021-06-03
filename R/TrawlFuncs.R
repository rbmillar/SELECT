
#****The core functions are ttfit() for alternate haul or trouser trawl data****
#****and ccfit() for covered codend data****************************************
#****and Rep.ttfit() for calculation of replicate estimate of overdispersion****
#*****from individual haul alternate or trouser trawl data**********************
#
#***********The sample program code in haddock.R, developed 9 July 2003*********
#*******************************************************************************
#*******************************************************************************
#The data are assumed to be in a matrix with 3 columns:
#Column 1:  Midpoint of lengthclass
#Column 3:  Numbers in control codend (trouser trawl) or cover
#Column 2:  Numbers in experimental codend
#
#EXAMPLES:  If the data matrix is called  catch  then:
#For a logistic fit to covered codend    ccfit(catch)
#For a Richards fit to covered codend    ccfit(catch,type="rich")
#For a logistic fit to trouser trawl     ttfit(catch)
#For a Richards fit to trouser trawl     ttfit(catch,type="rich")
#To fix split parameter to 0.5, say      ttfit(catch,psplit=0.5)
#
#Plots of fit and deviance resids are produced unless  plots=F  is specified.
#******************************************************************************


######################################
#Function for alternate haul analyses#
######################################
#' Fit logistic or Richards selection curve to paired hauls data.
#'
#' @description Fit logistic or Richards selection curve to paired hauls data. Pairing would typically be
#' via twin hauls, alternate hauls or trouser trawls
#' @param catch Matrix with 3 columns. These must be
#'  1) length (i.e., midpoint of length class);
#'  2) Freq in experimental gear;
#'  3) Freq in nonselective control gear
#' @param type Character. Use type="rich" to fit Richards curve
#' @param probs Numeric to specify retention probabilities for which corresponding
#'  lengths are required
#' @param psplit Numeric. If provided the split parameter is fixed at this value
#' @param plotlens Numeric. Lengths are which to predict retention probabilities.
#' @param details If TRUE, returns additional output
#' @param verbose If FALSE then messages are suppressed
#'
#' @return converged Convergence code from nlm minimizer
#' @return x Selection curve parameters, a, b, psplit and delta (if Richards curve)
#' @return l Log-likelihoods of fitted, null and full models
#' @return lens Lengths of 25%, 50% and 75% retention (by default)
#' @return sr Selection range
#' @return p psplit parameter
#' @return CF Overdispersion correction factors, using both deviance and Pearson chisq.
#'
#' @export

ttfit=function(catch=catchdat,type="logit",probs=c(0.25,0.5,0.75),
         psplit=NULL,x0=c(-10,0.4,0.6),delta=1.0,suff.big=3,nullfit=F,plots=T,
         plotfit=T,cex=0.8,mkh=0.07,error.bars=FALSE,plotlens=NULL,details=F,
         xlab="Length (cm)",ylab="Propn retained",verbose=T,delta.pen=0,
         main=c("Proportion of catch in large mesh codend","Deviance residuals")){
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]; ntotal=nfine+nwide
  nobs=length(ntotal[ntotal>0])
  if(!is.null(psplit)&verbose) cat("\n"," Fixed split, p= ",psplit)
  fullfithood=sum( dbinom(nwide,ntotal,ifelse(ntotal>0,nwide/ntotal,0),log=T) )
  nullfithood=sum( dbinom(nwide,ntotal,sum(nwide)/sum(ntotal),log=T) )

  if(length(plots)==1) plots=rep(plots,2)
  if(plots[1]) {
   propn=ifelse(ntotal>0.001,nwide/ntotal,2)
   xyticks=c(length(lenclass)-1,5,7)
   uniquelens=sort(unique(lenclass))
   AreLensUnique=(length(lenclass)==length(uniquelens))
   if(!error.bars)
    plot(lenclass[propn!=2], propn[propn!=2], pch = 5, mkh=mkh, las=1,
           type=ifelse(AreLensUnique,"b","p"),lab = xyticks, xlab = "", ylab = "",
           xlim = range(lenclass), ylim = c(0,1),cex=cex)
   else {
    lower.bnds=pmax(propn[propn!=2] + qnorm(0.1/2)*0.5/sqrt(ntotal)[propn!=2],0)
    upper.bnds=pmin(propn[propn!=2] - qnorm(0.1/2)*0.5/sqrt(ntotal)[propn!=2],1)
    errbar(lenclass[propn!=2],propn[propn!=2],lower.bnds,upper.bnds,incr=F,
          pch = 5,lab = xyticks,
          xlab = xlab,ylab = ylab,mkh=0.07,
          xlim = c(lenclass[1], lenclass[length(lenclass)]), ylim = c(0,1)) }
   title(xlab = xlab, ylab = ylab,main=main[1],cex=cex)}

  if(nullfit)
    {
    if(!is.null(psplit)) cselect=rep(psplit,length(lenclass))
      else cselect=rep(sum(nwide)/sum(ntotal),length(lenclass))
    }
  else
  {
  if(verbose)
    cat("\n","***NOTE: warning messages may occur as normal part of optimization***","\n")
  if(type=="logit")
   {
   if(!is.null(psplit)) {
    if(verbose) cat("\n","Fixed split p= ",psplit)
    Tfit=nlm(hood2par,x0[1:2],iterlim=1000,catch=catch,psplit=psplit);
    Pars=Tfit$est;
    Tcov=cov2par(Pars,catch,type="logit",p=psplit)
    select=lselect(Pars,lenclass)
    if(is.null(plotlens)) r=select else r=lselect(Pars,plotlens)
    cselect=psplit*select/(psplit*select + 1-psplit)
    Tlens=retentionlens(Pars,cov=Tcov$covar,probs=probs)
    p=cbind(psplit,NA) }
   else {
    #Tfit=optim(x0,hood3par,catch=catch); Pars=Tfit$par;
    #Tfit=nlm(hood3par,x0,iterlim=1000,catch=catch); Pars=Tfit$est;
    #Fit psplit on logit scale
    Tfit=nlm(hood3parlogit,c(x0[1:2],logit(x0[3])),iterlim=1000,catch=catch);
    Pars=c(Tfit$est[1:2],plogis(Tfit$est[3]));
    Tcov=cov3par(Pars,catch,type="logit")
    select=lselect(Pars,lenclass)
    if(is.null(plotlens)) r=select else r=lselect(Pars,plotlens)
    cselect=Pars[3]*select/(Pars[3]*select + 1-Pars[3])
    Tlens=retentionlens(Pars,cov=Tcov$covar[1:2,1:2],probs=probs)
    p=cbind(Pars[3],sqrt(Tcov$covar[3,3])) }
   }
  if(type=="rich")
   {
   if(!is.null(psplit)) {
    if(verbose) cat("\n","Fixed split p= ",psplit)
    Richhood2par=function(x,catch,psplit) richhood2par(x,catch,psplit)+delta.pen*log(x[3])^2
    Tfit=nlm(Richhood2par,c(x0[1:2],delta),iterlim=1000,catch=catch,psplit=psplit)
    Pars=Tfit$est
    if(verbose) cat("\n"," Likelihood of fitted model is ",
                    format(-richhood2par(Pars,catch,psplit)),"\n")
    select=lselect(Pars[1:2],lenclass)^(1/Pars[3])
    if(is.null(plotlens)) r=select else r=lselect(Pars[1:2],plotlens)^(1/Pars[3])
    cselect=psplit*select/(psplit*select + 1-psplit)
    Tcov=covrich(Pars,catch,npars=2,p=psplit)
    Tlens=retentionlens(Pars,cov=Tcov$covar[1:3,1:3],type="rich",probs=probs)
    p=cbind(psplit,NA)  }
   else {
    Richhood=function(x,catch) richhood(x,catch)+delta.pen*log(x[3])^2
    Tfit=nlm(Richhood,c(x0[1:2],delta,x0[3]),iterlim=1000,catch=catch); Pars=Tfit$est
    if(verbose) cat("\n"," Likelihood of fitted model is ",
                    format(-richhood(Pars,catch)),"\n")
    select=lselect(Pars[1:2],lenclass)^(1/Pars[3])
    if(is.null(plotlens)) r=select else r=lselect(Pars[1:2],plotlens)^(1/Pars[3])
    cselect=Pars[4]*select/(1-Pars[4] + Pars[4]*select)
    Tcov=covrich(Pars,catch)
    Tlens=retentionlens(Pars,cov=Tcov$covar[1:3,1:3],type="rich",probs=probs)
    p=cbind(Pars[4],sqrt(Tcov$covar[4,4])) }
     }
   #modelhood=-Tfit$value
   modelhood=-Tfit$min
   lhoods=cbind(c(modelhood,nullfithood,fullfithood),nobs-c(length(Pars),1,nobs))
   }
  yhat=ntotal*cselect
  Tdevres=devres(nwide,yhat,ntotal,suff.big)
  suff.dof=sum(Tdevres$suff.dat)-length(Pars)
  DevCF=sum(Tdevres$devres[Tdevres$suff.dat]^2)/suff.dof
  PearCF=sum(Tdevres$Pearson[Tdevres$suff.dat]^2)/suff.dof

  if(plots[1]&plotfit)
   lines(lenclass[order(lenclass)],cselect[order(lenclass)],type="l",lty=2)
  if(plots[2]) {
   plot(lenclass,Tdevres$devres,type=ifelse(AreLensUnique,"h","p"),las=1,
          lab=xyticks,xlab="",ylab="",cex=cex)
   abline(h=0)
   title(xlab=xlab,ylab="Deviance residual",main=main[2],cex=cex) }

  if(nullfit)
   {
    if(!is.null(psplit)) list(p=psplit)
      else list(p=sum(nwide)/sum(ntotal))
    }
  else
  if(details)
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr,p=p,
       CF=c(DevCF=DevCF,PearCF=PearCF),xcovar=Tcov$covar,lensr.covar=Tlens$covar,phi=cselect,
       r=r,devres=Tdevres$devres,suff.dat=Tdevres$suff.dat)
  else
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr,p=p,
        CF=c(DevCF=DevCF,PearCF=PearCF))}

######################################
#Function for covered codend analyses#
######################################
#' Fit logistic or Richards selection curve to covered-codend data.
#'
#' @param catch Matrix with 3 columns. These must be
#'  1) length (i.e., midpoint of length class);
#'  2) Freq in experimental gear;
#'  3) Freq in nonselective cover
#' @param type Character. Use type="rich" to fit Richards curve
#' @param probs Numeric to specify retention probabilities for which corresponding
#'  lengths are required
#' @param plotlens Numeric. Lengths are which to predict retention probabilities.
#' @param details If TRUE, returns additional output
#'
#' @return converged Convergence code from nlm minimizer
#' @return x Selection curve parameters, a, b, and delta (if Richards curve)
#' @return l Log-likelihoods of fitted, null and full models
#' @return lens Lengths of 25%, 50% and 75% retention (by default)
#' @return sr Selection range
#' @return CF Over-dispersion correction factors, using both deviance and Pearson chisq.
#'
#' @export
ccfit=function(catch=catchdat,type="logit",probs=c(0.25,0.5,0.75),x0=c(-10,0.3),
               delta=1.0,plots=T,suff.big=3,error.bars=F,plotlens=NULL,details=F,
               main=c("Proportion of catch in large mesh codend","Deviance residuals"),
               xlab="Length (cm)",ylab="Propn retained") {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]; ntotal=nfine+nwide
  nobs=length(ntotal[ntotal>0])
  fullfithood=sum( log(dbinom(nwide,ntotal,ifelse(ntotal>0,nwide/ntotal,0))) )
  nullfithood=sum( log(dbinom(nwide,ntotal,sum(nwide)/sum(ntotal))) )
  #cat("\n"," Log-likelihood of full model is ",format(fullfithood))
  #cat("\n"," Log-likelihood of null model is ",format(nullfithood))
  if(is.null(plotlens)) xlimits=range(lenclass)
    else xlimits=range(plotlens)

  if(length(plots)==1) plots=rep(plots,2)
  if(plots[1])
   {
   propn=ifelse(ntotal>0.001,nwide/ntotal,2)
   xyticks=c(length(lenclass)-1,10,7)
   uniquelens=sort(unique(lenclass))
   AreLensUnique=(length(lenclass)==length(uniquelens))
   if(!error.bars) {
    plot(lenclass[propn!=2], propn[propn!=2],pch = 5,lab = xyticks,
          type=ifelse(AreLensUnique,"b","p"),
          xlab = xlab, ylab = ylab,mkh=0.07,las=1,
          xlim = xlimits,ylim = c(0,1),main=main[1]) }
   else {
    lower.bnds=pmax(propn[propn!=2] + qnorm(0.05)*0.5/sqrt(ntotal)[propn!=2],0)
    upper.bnds=pmin(propn[propn!=2] - qnorm(0.05)*0.5/sqrt(ntotal)[propn!=2],1)
    error.bar(lenclass[propn!=2],propn[propn!=2],lower.bnds,upper.bnds,incr=F,
          pch = 5,lab = xyticks,
          xlab = xlab,ylab = ylab,mkh=0.07,
          xlim = xlimits, ylim = c(0,1), main=main[1]) }
   }
  if(type=="logit")
   {
    Tfit=nlm(cchood,x0,iterlim=100,catch=catch); Pars=Tfit$est;
    Tcov=cccov(Pars,catch)
    Tlens=retentionlens(Pars,Tcov$covar)
    select=lselect(Pars,lenclass)
    if(is.null(plotlens)) r=select else r=lselect(Pars,plotlens)
   }
  if(type=="rich")
   {
    Tfit=nlm(ccrichhood,c(x0,delta),iterlim=200,catch=catch); Pars=Tfit$est;
    Tcov=cccovrich(Pars,catch)
    Tlens=retentionlens(Pars,Tcov$covar,type="rich")
    select=lselect(Pars[1:2],lenclass)^(1/Pars[3])
    if(is.null(plotlens)) r=select else r=lselect(Pars[1:2],plotlens)^(1/Pars[3])
   }
  Tdevres=devres(nwide,ntotal*select,ntotal,suff.big)
  suff.dof=sum(Tdevres$suff.dat)-length(Pars)
  DevCF=sum(Tdevres$devres[Tdevres$suff.dat]^2)/suff.dof
  PCF=sum(Tdevres$Pearson[Tdevres$suff.dat]^2)/suff.dof

  modelhood=-Tfit$min
  lhoods=cbind(c(modelhood,nullfithood,fullfithood),nobs-c(length(Pars),1,nobs))
  if(plots[1]) {
    lines(lenclass[order(lenclass)],select[order(lenclass)],type="l",lty=2)
    abline(h=c(0.25,0.5,0.75),lty=3) }
  if(plots[2]) {
    plot(lenclass,Tdevres$devres,type=ifelse(AreLensUnique,"h","p"),main=main[2],
          xlim=xlimits,xlab=xlab,lab=xyticks,las=1,ylab="Deviance residual")
    abline(h=0) }
  ModelName=ifelse(type=="logit","logistic","Richards")
  cat("\n"," Log-likelihood of",ModelName,"model is ",modelhood,"\n")

  if(details)
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr,
     CF=c("DevCF"=DevCF,"PCF"=PCF), xcovar=Tcov$covar,lensr.covar=Tlens$covar,r=r,
     devres=Tdevres$devres,suff.dat=Tdevres$suff.dat)
  else
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr,
        CF=c("DevCF"=DevCF,"PCF"=PCF)) }

######################################################################
#Function for REP calculation for individual hauls alternate haul data
######################################################################
#'
#' Replicate estimation of over-dispersion for paired hauls data.
#'
#' @description  Replicate estimation of over-dispersion for logistic or Richards
#' selection curve fitted to stacked (by individual haul) paired hauls data.
#' Pairing would typically be via twin hauls, alternate hauls or trouser trawls.
#' This function applies the REP correction of over-dispersion in
#' Millar et al. (2004, Modelling between-haul variability in the size selectivity of trawls,
#'                Fisheries Research. 67: 171-181.)
#'
#' @param catch Matrix with 4 columns. These must be
#'  1) length (i.e., midpoint of length class);
#'  2) Freq in experimental gear; 3) Freq in nonselective control gear
#'  4) Haul ID
#' @param type Character. Use type="rich" to fit Richards curve
#' @param details If TRUE, returns additional output
#' @param verbose If FALSE then messages are suppressed
#'
#' @export
#'

#Replicate tows fits with individually estimated p's (if ind.psplits=T)
#Also fit null model (no selectivity)
#This version assumes that haul id is in column 4 of catchdat
Rep.ttfit=function(catch=catchdat,type="logit",x0=c(-9,0.3,0.5),
          delta=1.0,ind.psplit=T,suff.big=3,plots=F,details=T,verbose=T) {
  if(type!="logit" & type!= "rich") {
    cat("\nError: selection curve type must be logit or rich")
    return() }
  nobs=nrow(catch)
  towID=catch[,4]
  #if(sum(towID[-nobs]<=towID[-1])!=(nobs-1)) {
  #  cat("\nError: towID's in column 4 of catch data must be sorted")
  #  return() }
  uniqueTowID=unique(towID)
  ntows=length(unique(towID))
  if(verbose) cat("\nNumber of tows=",ntows)

  if(verbose) cat("\n","COMBINED TOWS FIT:")
  nobs=nrow(catch)
  nwide=catch[,2]
  nfine=catch[,3]
  lenorder=order(catch[,1]) #Use below to get plot right
  combfit=ttfit(catch[lenorder,],x0=x0,delta=delta,type=type,
                suff.big=suff.big,plots=plots,verbose=verbose)
  x=combfit$x
  if(verbose) {
    cat("\n Combined haul parameters:")
    if(type=="logit") cat("\n a, b and p:",x)
    if(type=="rich") cat("\n a, b, delta, and p:",x)
    if(ind.psplit)
      cat("\n","FIT OF COMBINED-HAULS CURVE TO REPLICATE TOWS WITH INDIVIDUAL SPLITS")
    else cat("\n","FIT OF COMBINED-HAULS CURVE TO REPLICATE TOWS WITH COMMON SPLIT")
  }

  suff.dat=rep(T,nobs)
  psplits=rep(0,ntows)

  llhood=c(0,0); names(llhood)=c("Fitted","Null") #log-likelihood for suff.big data
  llfull=c(0,0); names(llfull)=c("Fitted","Null") #log-likelihood for all data
  Dev=c(0,0) #Deviance for all data
  for(nonsel in c(F,T)){
  if(nonsel==F) k=1 else k=2;
  if(nonsel & verbose) {
    if(ind.psplit)
      cat("\n","NONSELECTIVE FIT WITH INDIVIDUAL TOW SPLITS")
    else cat("\n","NONSELECTIVE FIT USING A COMMON SPLIT")
    }

  P.suff=0; D.suff=0; LensUsed=0; TowsUsed=0
  #Loop over tow
  for(tow in 1:ntows) {
    indices=(towID==uniqueTowID[tow])
    lenclass=catch[indices,1]
    if(!nonsel) {
      if(type=="logit") r=lselect(x[1:2],lenclass)
      if(type=="rich") r=lselect(x[1:2],lenclass)^(1/x[3]) }
    else {rfit=r; r=1}
    nwide=catch[indices,2]
    nfine=catch[indices,3]
    n=nwide+nfine; p=ifelse(n>0,nwide/n,0); y=nwide
  #Calculate p and phi
    if(sum(n)>0) {
       #if(ind.psplit) psplit=sum(nwide)/sum(nwide+r*nfine) #Crude approx
       if(ind.psplit) {
          #psplit=nlm(hoodpsplit,0.5,catch=cbind(lenclass,nwide,nfine),r=r)$est
          logitpsplit=optim(0,hoodpsplit,catch=cbind(lenclass,nwide,nfine),r=r,
                       method="Brent",lower=-9,upper=9)$par
          psplit=plogis(logitpsplit) }
        else if(nonsel) psplit=sum(catch[,2])/sum(catch[,2]+catch[,3])
          else psplit=x[3]
      phi=psplit*r/((1-psplit)+psplit*r)
      yhat=n*phi
      PearsonStat=(y-yhat)^2/(n*phi*(1-phi))
      l=y*log(phi)+(n-y)*log(1-phi)
      lratio=ifelse(p>0,y*(log(p)-log(phi)),0)+
        ifelse(p<1,(n-y)*(log(1-p)-log(1-phi)),0)
      #Suff.dat based on fitted model
      if(!nonsel) suff.dat[indices]=(yhat>suff.big & (n-yhat)>suff.big)
      suff=suff.dat[indices]
      llfull[k]=llfull[k]+sum(l)
      Dev[k]=Dev[k]+2*sum(lratio)
      if(sum(suff)>0) {
        TowsUsed=TowsUsed+1
        llhood[k]=llhood[k]+sum(l[suff])
        P.suff=P.suff+sum(PearsonStat[suff])
        D.suff=D.suff+2*sum(lratio[suff])
        LensUsed=LensUsed+sum(suff) }
    }
    if(!nonsel) psplits[tow]=psplit
  }
  if(verbose) {
    cat("\n\n Using suff=",suff.big,": #lens=",LensUsed,"#tows=",TowsUsed,"and")
    cat("\n  Pearson Chisq=",round(P.suff,4),", Model dev=",round(D.suff,4),
              ", l=",round(llhood[k],4),"\n\n")
  }
  if(ind.psplit) dof=LensUsed-TowsUsed-2
        else if(nonsel) dof=LensUsed-1
          else dof=LensUsed-3
  if(!nonsel) {
    DevCF=D.suff/dof
    PearCF=P.suff/dof
    CorrectedMLEs=matrix(NA,5,2);
    CorrectedMLEs=rbind(combfit$lens,combfit$sr,combfit$p)
    CorrectedMLEs[,2]=CorrectedMLEs[,2]*sqrt(DevCF)
    dimnames(CorrectedMLEs)=list(c("L25","L50","L75","SR","p"),c("MLE","se"))
    if(verbose) {
      cat(" Using deviance, REP factor is ",
            round(D.suff,4),"/",dof,"=",round(DevCF,4),
          ", P-value=",1-pchisq(D.suff,dof),"\n")
      cat(" With correction applied: ","\n")
      print(CorrectedMLEs) }
  }
}
if(details)
  list(x=x,DevCF=DevCF,PearCF=PearCF,Dev=Dev[1],psplits=psplits,lsuff=llhood,l=llfull,
       Pars=CorrectedMLEs,r=rfit)
else
  list(DevCF=DevCF,PearCF=PearCF)
}

#'
#' Logistic selection curve
#'
#' @description The logistic selection curve function
#'
#' @param x Vector containing parameters a (intercept) and b (shape)
#' @param lenclass Vector of lengths
#' @returns Vector of retention probabilities
#'
#' @export
#'

lselect=function(x,lenclass) {
  expo=exp(pmin(500,x[1]+x[2]*lenclass))
  expo/(1+expo)
}

#'
#' Richards selection curve
#'
#' @description The Richards selection curve function
#'
#' @param x Vector containing parameters a (intercept), b (shape) and delta
#' @param lenclass Vector of lengths
#' @returns Vector of retention probabilities
#'
#' @export
#'
rselect=function(x,lenclass) {
  lselect(x[1:2],lenclass)^(1/x[3])
}


######################################
#Quadratic spline exploratory analysis#
######################################
#require(gamm4)
EDAfit=function(catch=catchdata,k=-1,fx=FALSE,m=c(2,1),cex=0.8) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3];
  ntotal=nfine+nwide; y=nwide/ntotal
  lenseq=seq(min(lenclass),max(lenclass),length=200)
  GAMfit=gam(y~s(lenclass,k=k,fx=fx,bs="bs",m=m),family=binomial,weights=ntotal)
  GAMpred=predict(GAMfit,data.frame(lenclass=lenseq),type="response")
  #Plot fit
  uniquelens=sort(unique(lenclass))
  AreLensUnique=(length(lenclass)==length(uniquelens))
  plot(lenclass, y, pch = 5,type=ifelse(AreLensUnique,"b","p"),
       xlab = "", ylab = "",xlim = range(lenclass), ylim = c(0,1),cex=cex)
  points(lenseq,GAMpred,type="l")
  list(GAMfit,GAMpred)
}

###############################################################################
##########END OF USER DESIGNED FUNCTIONS#######################################
###############################################################################

#******************************************************************************
#These next functions are called by ttfit(), ccfit(), and Rep.ttfit()
#and would not normally be called by the user.
#******************************************************************************

logit=function(p) log(p/(1-p))

lselect=function(x,lenclass) {
  expo=exp(pmin(500,x[1]+x[2]*lenclass))
  expo/(1+expo) }

rselect=function(x,lenclass) {
  lselect(x[1:2],lenclass)^(1/x[3])
}

cchood=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  select=lselect(x,lenclass)
  #-sum(nwide*log(select) + nfine*log(1-select))
  -sum(dbinom(nwide,nwide+nfine,select,log=T))
  }

ccrichhood=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  select=lselect(x[1:2],lenclass)^(1/x[3])
  #-sum( nwide*ifelse(select>0,log(select),-1e+06) +
  #    nfine*ifelse(select<1,log(1-select),-1e+06) )
  -sum(dbinom(nwide,nwide+nfine,select,log=T))
  }

hood2par=function(x,catch,psplit) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  expo=exp(x[1]+x[2]*lenclass)
  cselect=psplit*expo/( 1-psplit + expo)
  -sum(dbinom(nwide,nwide+nfine,cselect,log=T)) }

hood3par=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  expo=exp(pmin(500,x[1]+x[2]*lenclass))
  cselect=x[3]*expo/( (1-x[3]) + expo)
  -sum(dbinom(nwide,nwide+nfine,cselect,log=T)) }

hood3parlogit=function(x,catch) {
  psplit=plogis(x[3])
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  expo=exp(pmin(500,x[1]+x[2]*lenclass))
  cselect=psplit*expo/( (1-psplit) + expo)
  -sum(dbinom(nwide,nwide+nfine,cselect,log=T)) }

hoodpsplit=function(logitpsplit,catch,r,r2=1) {
  psplit=plogis(logitpsplit)
  nwide=catch[,2]; nfine=catch[,3]
  cselect=psplit*r/((1-psplit)*r2+psplit*r)
  -sum( nwide*ifelse(cselect>0,log(cselect),-500) +
      nfine*ifelse(cselect<1,log(1-cselect),-500) ) }

richhood=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  #select=lselect(x[1:2],lenclass)^(1/x[3])
  select=rselect(x,lenclass)
  cselect=x[4]*select/(1-x[4] + x[4]*select)
  -sum(dbinom(nwide,nwide+nfine,cselect,log=T)) }

richhood2par=function(x,catch,psplit) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  #select=lselect(x[1:2],lenclass)^(1/x[3])
  select=rselect(x,lenclass)
  cselect=psplit*select/(1-psplit + psplit*select)
  -sum(dbinom(nwide,nwide+nfine,cselect,log=T)) }

#Returns the Pearson and deviance residuals
devres=function(y,yhat,n,suff.big=3,verbose=F) {
  if( any((n==0)&(y>0 | yhat>0)) ) stop("Wrong data in function devres")
  if( any((yhat==0 & y >0) | (yhat==n & y<n)) )
    stop("Impossibility in function devres")
  suff.dat=(yhat>suff.big & (n-yhat)>suff.big)
  p=ifelse(n>0,y/n,0.5); phat=ifelse(n>0,yhat/n,0.5)
  sign=ifelse(y>=yhat,1,-1)
  Pearson=(y-yhat)/ifelse(n*phat*(1-phat)>0,sqrt(n*phat*(1-phat)),1)
  l=ifelse(y>0,y*(log(p)-log(phat)),0) +
     ifelse(y<n,(n-y)*(log(1-p)-log(1-phat)),0)
  if(verbose) {
    cat("\n"," Pearson Chisq=",round(sum(Pearson^2),4),
        ", Dev=",round(2*sum(l),4)," #lens used=",sum(n>0),sep="")
    cat("\n"," Pearson Chisq=",round(sum(Pearson[suff.dat]^2),4),
        ", Dev=",round(2*sum(l[suff.dat]),4),
        " #lens used=",sum(suff.dat)," (Expected count >",suff.big,")","\n",sep="") }
  list(Pearson=Pearson,devres=sign*sqrt(2*l),suff.dat=suff.dat) }

retentionlens=function(x,covar,type="logit",probs=c(0.25,0.5,0.75)) {
  np=length(probs)
  if(type=="logit") {
    rlens=( log(probs/(1-probs)) -x[1] ) / x[2]
    srange=2*log(3)/x[2]
    if(!missing(covar)) {
      derivs=matrix(0,nrow=2,ncol=np+1)
      for(i in 1:np) derivs[,i]=c(-1/x[2],-rlens[i]/x[2])
      derivs[,np+1]=c(0,-srange/x[2])
      rlencovar= t(derivs) %*% covar %*% derivs
      lens=matrix(c(rlens,sqrt(diag(rlencovar)[1:np])),nrow=np,byrow=F)
      sr=c(srange,sqrt(rlencovar[np+1,np+1]))
      return(list(lens=lens,sr=sr,covar=rlencovar)) } }
  if(type=="rich") {
    work=(probs^x[3])/(1-probs^x[3]); rlens=(log(work)-x[1])/x[2]
    worksr=(c(0.25,0.75)^x[3])/(1-c(0.25,0.75)^x[3])
    srange=(log(worksr[2])-log(worksr[1]))/x[2]
    if(!missing(covar)) {
     derivs=matrix(0,nrow=3,ncol=np+1)
     for(i in 1:np)
      derivs[,i]=c(-1/x[2],-rlens[i]/x[2],log(probs[i])/(x[2]*(1-probs[i]^x[3])))
      derivs[,np+1]=c(0,-srange/x[2],
                     (log(0.75)/(1-0.75^x[3])-log(0.25)/(1-0.25^x[3]))/x[2])
     rlencovar= t(derivs) %*% covar %*% derivs
     lens=matrix(c(rlens,sqrt(diag(rlencovar)[1:np])),nrow=np,byrow=F)
     rownames(lens)=c("L25","L50","L75"); colnames(lens)=c("MLE","se")
     sr=c(srange,sqrt(rlencovar[np+1,np+1])); names(sr)=c("MLE","se")
     return(list(lens=lens,sr=sr,covar=rlencovar)) } }
   }

#**************************************************************************
#**************************************************************************
#For historical reasons, the covariance matrix of the ML estimates are
#obtained from expected Fisher information, using the functions below
#**************************************************************************
#**************************************************************************

#***If problems occur with trouser trawl cov functions then numerical******
#***accuracy enhancements used in cccovrich will need to be added**********

#Information matrix and diagnostic calculations for fixed split fit.
cov2par=function(x,catch,type="logit",p=0.5) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  if(type=="logit") {
    mu=p*exp(eta)/(1-p+exp(eta))
    dmudeta=p*(1-p)*exp(eta)/( (1-p+exp(eta))^2 ) }
  if(type=="cl") {
    dblexpo=exp(-exp(eta))
    mu=1-(1-p)/(1-p*dblexpo)
    dmudeta=(1-p)*p*exp(eta)*dblexpo/ ( (1-p*dblexpo)^2 ) }
  if(type=="logit" | type=="cl") {
    info=matrix(0,nrow=2,ncol=2)
    info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
    info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
    info[2,1]=info[1,2]
    info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
    covar=Solve(info)
    list(covar=covar)  }
  else return(NA)  }

#Information matrix and diagnostic calculations for estimated p fit.
cov3par=function(x,catch,type="logit") {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  if(type=="logit") {
    mu=(x[3]*exp(eta))/(1-x[3]+exp(eta))
    dmudeta= x[3]*(1-x[3])*exp(eta) / ( (1-x[3]+exp(eta))^2 )
    dmudp= exp(eta)*(1+exp(eta)) / ( (1-x[3]+exp(eta))^2 ) }
  if(type=="cl") {
    dblexpo=exp(-exp(eta))
    mu=1-(1-x[3])/(1-x[3]*dblexpo)
    dmudeta=(1-x[3])*x[3]*exp(eta)*dblexpo/ ( (1-x[3]*dblexpo)^2 )
    dmudp=(1-dblexpo)/( (1-x[3]*dblexpo)^2 ) }
  if(type=="logit" | type=="cl") {
    info=matrix(0,nrow=3,ncol=3)
    info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
    info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
    info[1,3]=sum( ntotal*dmudp*dmudeta / (mu*(1-mu)) )
    info[2,1]=info[1,2]
    info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
    info[2,3]=sum( ntotal*lenclass*dmudeta*dmudp / (mu*(1-mu)) )
    info[3,1]=info[1,3]
    info[3,2]=info[2,3]
    info[3,3]=sum( ntotal*dmudp^2 / (mu*(1-mu)) )
    covar=Solve(info,tol=1e-12)
    list(covar=covar)  }
  else return(NA)  }


#Covariance matrix for logistic fit to covered codend data
cccov=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  mu=exp(eta)/(1+exp(eta))
  dmudeta=mu/(1+exp(eta))
  info=matrix(0,nrow=2,ncol=2)
  info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
  info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
  info[2,1]=info[1,2]
  info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
  covar=Solve(info)
  list(covar=covar) }

#Covariance matrix for Richard's fit to covered codend data
cccovrich=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  nu=exp(eta)/(1+exp(eta)); mu=nu^(1/x[3])
  eta <- ifelse(mu == 1, Inf, eta) # To avoid numerical problems below
  dmudeta=(1/x[3])*mu/(1+exp(eta))
  dmudgamma=-x[3]^(-2)*mu*log(nu)
  info=matrix(0,nrow=3,ncol=3)
  info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) ,na.rm=T)
  info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) ,na.rm=T)
  info[1,3]=sum( ntotal*dmudgamma*dmudeta / (mu*(1-mu)) ,na.rm=T)
  info[2,1]=info[1,2]
  info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) ,na.rm=T)
  info[2,3]=sum( ntotal*lenclass*dmudeta*dmudgamma / (mu*(1-mu)) ,na.rm=T)
  info[3,1]=info[1,3]
  info[3,2]=info[2,3]
  info[3,3]=sum( ntotal*dmudgamma^2 / (mu*(1-mu)) ,na.rm=T)
  covar=Solve(info)
  list(covar=covar) }

#Covariance matrix for Richard's fit to trouser trawl data
covrich=function(x,catch,npars=3,p=0.5) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  if(npars==3) p=x[4]
  nu=exp(eta)/(1+exp(eta)); r=nu^(1/x[3]); mu=p*r/(1-p+p*r)
  dmudr=p*(1-p)/((1-p+p*r)^2)
  drdeta=(1/x[3])*r/(1+exp(eta)); dmudeta=dmudr*drdeta
  drdgamma=-x[3]^(-2)*r*log(nu);  dmudgamma=dmudr*drdgamma
  dmudp=r/((1-p+p*r)^2)
  info=matrix(0,nrow=4,ncol=4)
  info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
  info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
  info[1,3]=sum( ntotal*dmudgamma*dmudeta / (mu*(1-mu)) )
  info[1,4]=sum( ntotal*dmudeta*dmudp / (mu*(1-mu)) )
  info[2,1]=info[1,2]
  info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
  info[2,3]=sum( ntotal*lenclass*dmudeta*dmudgamma / (mu*(1-mu)) )
  info[2,4]=sum( ntotal*lenclass*dmudeta*dmudp / (mu*(1-mu)) )
  info[3,1]=info[1,3]
  info[3,2]=info[2,3]
  info[3,3]=sum( ntotal*dmudgamma^2 / (mu*(1-mu)) )
  info[3,4]=sum( ntotal*dmudgamma*dmudp / (mu*(1-mu)) )
  info[4,1]=info[1,4]
  info[4,2]=info[2,4]
  info[4,3]=info[3,4]
  info[4,4]=sum( ntotal*dmudp*dmudp / (mu*(1-mu)) )
  if(npars==3) covar=Solve(info) else covar=Solve(info[1:3,1:3])
  list(covar=covar) }

#Robust version of solve()
Solve=function(M,...) {
  inverse=try(solve(M))
  if(class(inverse)=="try-error") {
    #eps=1
    #cat("\n***Covariances are invalid***\n")
    #inverse=solve(M+diag(eps,nrow(M)),...)
    #inverse=diag(1,nrow(M))
    inverse=diag(NA,nrow(M)) }
  return(inverse)
}








