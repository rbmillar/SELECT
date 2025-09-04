## Model evaluation and helper functions
# ModelCheck, Estimates, PlotCurves
# nllhood, calcOD, StartValues, mchoose, deltamethod, Raw2Tots & SELECT.FORMAT

#===============================================================================
#' Provide a diagnostic summary of the SELECT model fit
#' @description Provides log-likelihoods, deviances,
#' over-dispersion correction factors (if minE>0), and a plot of deviance residuals.
#' @export
ModelCheck=function(fit,minE=0,xlab="Length (cm)",ylab = "Propn in exptl gear",
                    print.out=T,plots=c(T,T),plotlens=NULL,pex=1,...) {
  lens=fit$Data[,1]; nlens=length(lens)
  Meshsize=fit$Meshsize; nmeshes=length(Meshsize)
  O=fit$Data[,-1]; #Matrix of observed counts
  Const=mchoose(O,log=T) #Constant for log-likelihoods
  Ototals=apply(O,1,sum,na.rm=TRUE)
  Opropns=O/Ototals
  NullPropns=matrix(apply(O,2,sum,na.rm=TRUE),nrow=nlens,ncol=nmeshes,byrow=T)/sum(O)
  r=propncurves(fit$rtype) #Get propn catch curve function
  rmatrix=outer(lens,Meshsize,r,fit$par)
  rmatrix[is.na(O)]=NA #No fitted retention for missing meshsizes
  rmatrix=t(t(rmatrix)*fit$rel.power)
  phi=rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
  E=apply(O,1,sum,na.rm=TRUE)*phi
  Pearson.resids=(O-E)/sqrt(E)
  Pearson.chisq=sum(Pearson.resids^2,na.rm=TRUE)
  wk=O*log(O/E); wk[is.na(wk)]=0
  Dev.resids=sign(O-E)*sqrt(2*(E-O+wk))
  Deviance=sum(Dev.resids^2,na.rm=TRUE)
  #full.l=sum(-O+O*log(O),na.rm=TRUE)
  #null.l=sum(-null.E+O*log(null.E),na.rm=TRUE)
  #model.l=sum(-E+O*log(E),na.rm=TRUE)
  full.l=Const+sum(O*log(Opropns),na.rm=TRUE)
  null.l=Const+sum(O*log(NullPropns),na.rm=TRUE)
  model.l=Const+sum(O*log(phi),na.rm=TRUE)
  NonZeroDat=O[apply(O,1,sum,na.rm=TRUE)>0,]
  npar=length(fit$par)
  aic=-2*model.l+2*npar
  bic=-2*model.l+log(fit$nobs)*npar
  dof=nrow(NonZeroDat)*(nmeshes-1)-npar-sum(is.na(NonZeroDat))
  Deviance.CF=Deviance/dof; Pearson.CF=Pearson.chisq/dof
  out1a=cbind(null.l,model.l,full.l,npar,AIC=aic,BIC=bic)
  out1b=cbind(Deviance,Pearson.chisq,dof,Deviance.CF,Pearson.CF)
  outlist=list(stats=out1a,gof=out1b,fit=E)
  #If n cells for a given length have freq>=minE, it contributes max(0,n-1) dof
  if(minE>0) {
    Index=(E>minE)
    #But also need to exclude a cell if it is the only one in a row with freq>=minE
    RowIndex=(1:nrow(E))[apply(Index,1,sum,na.rm=TRUE)==1]
    Index[RowIndex,]=FALSE
    dof=sum(pmax(0,apply(Index,1,sum,na.rm=TRUE)-1))-npar
    Pearson.chisq=sum((Pearson.resids^2)[Index],na.rm=TRUE)
    Deviance=sum((Dev.resids^2)[Index],na.rm=TRUE)
    Deviance.CF=Deviance/dof
    Pearson.CF=Pearson.chisq/dof
    out2=cbind(Deviance,Pearson.chisq,dof,Deviance.CF,Pearson.CF)
  }
  #plots argument controls plotting of both deviance and fits

  if(length(plots)==1) plots=c(plots,plots)
  xyticks=c(nlens-1,10,7)
  AreLensUnique=(length(lens)==length(unique(lens)))
  if(plots[1]) { #Plot deviance residuals
  if(nmeshes>2&AreLensUnique) {
    plot(1,1,xlim=range(lens),xlab=xlab,ylab="Mesh size",
         ylim=range(Meshsize)+(pex/50)*c(-1,1)*(max(Meshsize)-min(Meshsize)),
         yaxt="n",type="n",...)
    axis(2,Meshsize,Meshsize)
    for(i in 1:nlens)
      for(j in 1:nmeshes)
        points(lens[i],Meshsize[j],pch=ifelse(Dev.resids[i,j]>0,16,1),
               cex=3*abs(Dev.resids[i,j])*pex/(abs(max(Dev.resids)))) }
  else
  if(nmeshes==2) {
    #Plot deviance residuals
    Dev.resids.len=sign(Dev.resids[,2])*sqrt(apply(Dev.resids^2,1,sum))
    plot(lens,Dev.resids.len,type=ifelse(AreLensUnique,"h","p"),
         xlab=xlab,ylab="Deviance residuals",lab=xyticks,...)
    abline(h=0) }
  } #End of plotting deviance residuals
  if(!is.null(plotlens)) {
    rmatrix=outer(plotlens,fit$Meshsize,r,fit$par)
    rmatrix=t(t(rmatrix)*fit$rel.power)
    phi=rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
  }
  if(is.null(plotlens)) plotlens=lens
  rownames(phi)=plotlens
  outlist=c(outlist,list(phi=phi))
  if( plots[2] & nmeshes==2){ #Plot the fit to propns
    Data=fit$Data; pwr=fit$rel.power
	  plot(lens,Data[,3]/(Data[,2]+Data[,3]),type=ifelse(AreLensUnique,"b","p"),
             ylim=c(0,1),xlab=xlab,ylab=ylab,lab=xyticks,...)
	  lines(plotlens,phi[,2],lty=2,...)
	}
  if(print.out) {
    cat("Model fit:\n"); print(out1a[1,])
    cat("GOF:\n"); print(out1b[1,])
    if(minE>0) {
      cat("\nCorrection factors from cells with expected count >",minE,":\n")
      print(out2[1,])
	    outlist=c(outlist,CF=out2) } }
  invisible(outlist)
}


#===============================================================================
#' Returns selectivity parameters
#' @description Calculates the estimates and standard errors for
#' parametric selectivity curves fitted by SELECT.
#' The user can edit this function to add new
#' calculations for any user-supplied retention curves.
#' @export
Estimates=function(fit,OD=NULL) {
  if(!is.null(OD)&!is.numeric(OD)) stop("\n OD must be numeric")
  x=fit$par;
  if(is.null(fit$hess))
    stop("\n SELECT fit does not include a hessian. Perhaps you used fit=F.\n")
  varx=solve(fit$hess)
  names=c("Mode(mesh1)","Std dev.(mesh1)")
  rtype=fit$rtype
  switch( #"dc" and "cc" estimates have identical parameter sets
    {rtype=ifelse(rtype=="cc.logistic","dc.logistic",rtype)
     rtype=ifelse(rtype=="cc.richards","dc.richards",rtype)},
    "dc.norm.loc"={ pars=x; varpars=varx },
    "dc.norm.sca"=, "dc.normal"={ pars=x; varpars=varx },
    "dc.gamma"={
      pars=c((x[1]-1)*x[2],sqrt(x[1]*x[2]^2))
      varpars=deltamethod(list(~(x1-1)*x2,~x2*sqrt(x1)),x,varx,ses=F)},

    "dc.lognorm"={
      pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)))
      varpars=deltamethod(list(~exp(x1-x2^2),
               ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),x,varx,ses=F)},
    "dc.binorm.sca"=, "dc.binormal"={
      pars=c(x[1:4],exp(x[5])/(1+exp(x[5])))
      names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                    "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
      varpars=deltamethod(list(~x1,~x2,~x3,~x4,~exp(x5)/(1+exp(x5))),
                                                x,varx,ses=F)},
    "dc.bilognorm"={
      pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)),
             exp(x[3]-x[4]^2),sqrt(exp(2*x[3]+x[4]^2)*(exp(x[4]^2)-1)),
             exp(x[5])/(1+exp(x[5])))
      names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                    "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
      varpars=deltamethod(
        list(~exp(x1-x2^2),~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1)),
             ~exp(x3-x4^2),~sqrt(exp(2*x3+x4^2)*(exp(x4^2)-1)),
             ~exp(x5)/(1+exp(x5))),x,varx,ses=F)},
	"ec.logistic"={
      pars=c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
      names=c("L50","SR","p")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
                          x,varx,ses=F)},
	"dc.logistic"={
      pars=c(-x[1]/x[2],2*(log(3))/x[2])
      names=c("L50","SR")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2),x,varx,ses=F)},
    "ec.richards"={
      delta=exp(x[3])
      pars=c((log(0.5^delta/(1-0.5^delta))-x[1])/x[2],
             (log(0.75^delta/(1-0.75^delta))-log(0.25^delta/(1-0.25^delta)))/x[2],
             delta,exp(x[4])/(1+exp(x[4])))
      names=c("L50","SR","delta","p")
      varpars=deltamethod(list(
        ~(log(0.5^exp(x3)/(1-0.5^exp(x3)))-x1)/x2,
        ~(log(0.75^exp(x3)/(1-0.75^exp(x3))))/x2
        -(log(0.25^exp(x3)/(1-0.25^exp(x3))))/x2,
        ~exp(x3),~exp(x4)/(1+exp(x4))),x,varx,ses=F)},
	"dc.richards"={
      delta=exp(x[3])
      pars=c((log(0.25^delta/(1-0.25^delta))-x[1])/x[2],
             (log(0.5^delta/(1-0.5^delta))-x[1])/x[2],
             (log(0.75^delta/(1-0.75^delta))-x[1])/x[2],
             (log(0.75^delta/(1-0.75^delta))-log(0.25^delta/(1-0.25^delta)))/x[2],
             delta)
      names=c("L25","L50","L75","SR","delta")
      varpars=deltamethod(list(
        ~(log(0.25^exp(x3)/(1-0.25^exp(x3)))-x1)/x2,
        ~(log(0.5^exp(x3)/(1-0.5^exp(x3)))-x1)/x2,
        ~(log(0.75^exp(x3)/(1-0.75^exp(x3)))-x1)/x2,
        ~(log(0.75^exp(x3)/(1-0.75^exp(x3))))/x2
        -(log(0.25^exp(x3)/(1-0.25^exp(x3))))/x2,~exp(x3)),x,varx,ses=F)},

    ###################Deprecated options#######################################
	"ec.logistic.L50SR"={
      pars=c(x[1],x[2],plogis(x[3]))
      names=c("L50","SR","p")
      varpars=deltamethod(list(~x1,~x2,~exp(x3)/(1+exp(x3))),x,varx,ses=F)},
    "ec.richards.L50SR"={
      pars=c(x[1],x[2],exp(x[3]),plogis(x[4]))
      print(pars)
      names=c("L50","SR","delta","p")
      varpars=deltamethod(list(~x1,~x2,~exp(x3),~exp(x4)/(1+exp(x4))),x,varx,ses=F)},
	"dc.logistic.L50SR"={ pars=x; varpars=varx; names=c("L50","SR") },
    "dc.logistic.ab"={ #If using a,b parameterization (old code)
      pars=c(-x[1]/x[2],2*(log(3))/x[2])
      names=c("L50","SR")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2),x,varx,ses=F)},
	"dc.richards.L50SR"={
      pars=c(x[1],x[2],exp(x[3]))
      names=c("L50","SR","delta")
      varpars=deltamethod(list(~x1,~x2,~exp(x3)),x,varx,ses=F)},
    "dc.richards.ab"={ #If using a,b parameterization (old code)
      delta=exp(x[3])
      pars=c((log(0.5^delta/(1-0.5^delta))-x[1])/x[2],
             (log(0.75^delta/(1-0.75^delta))-log(0.25^delta/(1-0.25^delta)))/x[2],
             delta)
      names=c("L50","SR","delta")
      varpars=deltamethod(list(
        ~(log(0.5^exp(x3)/(1-0.5^exp(x3)))-x1)/x2,
        ~(log(0.75^exp(x3)/(1-0.75^exp(x3))))/x2
        -(log(0.25^exp(x3)/(1-0.25^exp(x3))))/x2,~exp(x3)),x,varx,ses=F)},
    "ec.logistic.ab"={
      pars=c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
      names=c("L50","SR","p")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
                          x,varx,ses=F)},
    "ec.richards.ab"={
      delta=exp(x[3])
      pars=c((log(0.5^delta/(1-0.5^delta))-x[1])/x[2],
             (log(0.75^delta/(1-0.75^delta))-log(0.25^delta/(1-0.25^delta)))/x[2],
             delta,exp(x[4])/(1+exp(x[4])))
      names=c("L50","SR","delta","p")
      varpars=deltamethod(list(
        ~(log(0.5^exp(x3)/(1-0.5^exp(x3)))-x1)/x2,
        ~(log(0.75^exp(x3)/(1-0.75^exp(x3))))/x2
        -(log(0.25^exp(x3)/(1-0.25^exp(x3))))/x2,
        ~exp(x3),~exp(x4)/(1+exp(x4))),x,varx,ses=F)},
  stop(paste0('SELECT errror message: ',fit$rtype, ' not recognised.\n',
      'Possible direct comparison types are , "norm.loc", "normal", "gamma",
        "lognorm", "logistic", "richards, "binormal", and "bilognorm". \n',
      'Possbile covered-codend and alternative hauls types are
                   "logistic" and "richards" \n'))
  )#End of switch
  std.errors=sqrt(diag(varpars))
  estimates=cbind(pars,std.errors)
  colnames(estimates)=c("par","raw s.e.")
  if(!is.null(OD)) {
    estimates=cbind(estimates,sqrt(OD)*std.errors)
    colnames(estimates)[3]="adj s.e." }
  rownames(estimates)=names
  return(estimates) }



#===============================================================================
#' Plots the fitted curves and returns fitted values at specified lengths
#' @description Plot the fitted selectivity curve(s),
#' and return the fitted values if assigned to an object.
#' @param fit Fitted SELECT model object
#' @param plotlens Vector of lengths at which to calculate the fitted curves.
#' @param Meshsize Vector of mesh sizes if desired to calculate fitted values
#' at other meshsizes. Primarily for design type `dc`
#' @param rel.power Relative power for use with design type `dc`
#' @param standardize Standardize the fitted values to a maximum value of 1.
#' Default is F for design types `cc` and `ec`, and T for design type `dc`.
#' @param Master For design type `dc`, if TRUE (default) then the master curve
#' is used. This assumes geometric similarity and calculates retention as
#' a function of length/meshsize.
#' This should be set to F with `norm.loc` curves as these are not
#' geometrically similar.
#' @return A dataframe with fitted values at the specified lengths.
#' @export
PlotCurves=function(fit,plotlens=NULL,Meshsize=NULL,rel.power=NULL,
                    standardize=NULL,Master=T,plot.out=T,ylim=c(0,1),
                    xlab="Length (cm)",ylab="Retention curve",type="l",...) {
  s=selncurves(fit$rtype) #Get selection curve function
  lgthName=fit$var.names[1]
  dtype=substr(fit$rtype,1,2)
  if(is.null(plotlens)) plotlens=fit$Data[,1]
  if(is.null(Meshsize)) Meshsize=fit$Meshsize
  if(is.null(rel.power)) pwr=fit$rel.power
  if(is.null(standardize)) standardize=ifelse(dtype=="dc",T,F)
  smatrix=outer(plotlens,Meshsize,s,fit$par)
  smatrix=t(t(smatrix)*pwr)
  if(standardize) smatrix=smatrix/max(smatrix)
  Df=as.data.frame(cbind(plotlens,smatrix))
  colnames(Df)=fit$var.names
  if(dtype %in% c("cc","ec") & plot.out) {
    plot(plotlens,smatrix[,2],ylim=ylim,xlab=xlab,ylab=ylab,type=type,...)
    abline(h=c(0.25,0.5,0.75),lty=3) }
  if(dtype == "dc" & !Master & plot.out)
    matplot(plotlens,smatrix,ylim=ylim,xlab=xlab,ylab=ylab,type=type,...)
  if(dtype == "dc" & Master) {
    colnames(Df)=c(fit$var.names[1],Meshsize)
    Df=pivot_longer(Df,cols=!all_of(lgthName),values_to="r",names_to="Msize")
    Df=as.data.frame(Df)
    Df$Msize=as.numeric(Df$Msize)
    Df$v=Df[,lgthName]/Df$Msize
    Df=Df[order(Df$v),]
    if(plot.out) plot(r~v,type=type,xlab="Length/Meshsize",ylab=ylab,data=Df) }
  invisible(Df)
}



#===============================================================================
#' Negative log-likelihood of the SELECT model (excluding constants)
#' @description The general negative log-likelihood function. Provided for
#' completeness, but not intended for user use.
#' @export
nllhood=function(theta,Data,Meshsize,r,rel.power,penalty.func) {
  lens=Data[,1]; Counts=Data[,-1]
  rmatrix=outer(lens,Meshsize,r,theta)
  rmatrix[is.na(Counts)]=0 #No fitted retention for NA counts
  rmatrix=t(t(rmatrix)*rel.power)
  phi=rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
  nll=-sum(Counts[Counts>0]*log(phi[Counts>0]),na.rm=T)
  #nll=sum(byrow_dmultinom(Counts,prob=phi,log=T))
  nll=nll+penalty.func(theta)
  return(nll)
}

#===============================================================================
#' Calculate quasi-AIC
#' @description Adjusts the AIC of a fit for the given overdisperion `OD`.
#' AICc is used if `correct=T`.
#' @export
qAIC=function(fit,OD,correct=F) {
  Correction=0
  if(correct) {
    k=fit$rank
    Correction=2*k*(k+1)/(fit$df.resid-1)
  }
  qAIC=-2*logLik(fit)/OD+2*(attributes(logLik(fit))$df)+Correction
}


#===============================================================================
#' Function to calculate overdispersion
#' @description Returns OD estimate from matrices of observed and expected counts.
#' If the two matrices have equal row sum than the degrees of freedom is adjusted
#' for this constraint.
#' @param O Matrix of observed counts
#' @param E Matrix of expected (fitted) counts
#' @param minE Calculates OD only over cells with at least minE expected count.
#' @param verbose If FALSE then the message is suppressed.
#' @export
calcOD=function(O,E,npar,minE=1,verbose=T) {
  if(!identical(dim(O),dim(E))) stop("O and E must be the same dimension")
  EqualRowSums={ max(abs(apply(O,1,sum)-apply(E,1,sum)))<1e-6 }
  if(verbose & EqualRowSums)
    cat("\nNOTE: O and E have equal row sums, so dof will be reduced by
         the number of rows used") else
    cat("\nNOTE: O and E have unequal row sums, i.e., non-SELECT fit")
  dof.adj=ifelse(EqualRowSums,1,0)

  Pearson.resids=(O-E)/sqrt(E)
  wk=O*log(O/E); wk[is.na(wk)]=0
  Dev.resids=sign(O-E)*sqrt(2*(E-O+wk))
  #If n cells for a given length have freq>=minE, it contributes max(0,n-1) dof
  Index=(E>minE)
  #But also need to exclude a cell if it is the only one in a row with freq>=minE
  RowIndex=(1:nrow(E))[apply(Index,1,sum,na.rm=TRUE)==1]
  Index[RowIndex,]=FALSE
  dof=sum(pmax(0,apply(Index,1,sum,na.rm=TRUE)-dof.adj))-npar
  Pearson.chisq=sum((Pearson.resids^2)[Index],na.rm=TRUE)
  Deviance=sum((Dev.resids^2)[Index],na.rm=TRUE)
  Deviance.CF=Deviance/dof
  Pearson.CF=Pearson.chisq/dof
  ODstats=cbind(Deviance,Pearson.chisq,dof,Deviance.CF,Pearson.CF)
  return(ODstats)
  }

#===============================================================================
#' Initial parameter values for SELECT fits
#' @description Tf SELECT is not provided with explicit initial parameter values
#' values `x0` then crude data-driven starting values are generated.
#' The user will have to provide explicit `x0` if these default value do not
#' achieve convergence or a sensible fit.
#' @details The starting values are based on the median and standard deviation
#' of the lengths in the codend if the design type is `cc` or `ec`,
#' or of the smallest mesh size if design type is `dc`.
#'
#'
#' @export
StartValues=function(rtype,Data) {
  if(substr(rtype,1,2)=="cc" | substr(rtype,1,2)=="ec") {
    CodendMed=median(rep(Data[,1],Data[,3]))
    CodendSd=sd(rep(Data[,1],Data[,3]))
    b0=2*log(3)/CodendSd; a0=-b0*CodendMed
    switch(substr(rtype,1,6),
         "cc.log"={ c(a0,b0) },
         "cc.ric"={ c(a0,b0,0) },
         "ec.log"={ c(a0,b0,0) },
         "ec.ric"={ c(a0,b0,0,0) },
         stop("Please provide a value of x0 (initial parameter values")
    ) }
  else {
    Mesh1Med=median(rep(Data[,1],Data[,2]))
    Mesh1Var=var(rep(Data[,1],Data[,2]))
    Mesh1SD=sqrt(Mesh1Var)
    CV=Mesh1SD/Mesh1Med
    mu.sd=c(Mesh1Med,Mesh1SD)
    lnorm.mu.sd=c(log(Mesh1Med),CV)
    switch(substr(rtype,1,6),
           "dc.nor"={ mu.sd },
           "dc.log"={ lnorm.mu.sd },
           "dc.gam"={ c(Mesh1Med^2/Mesh1Var,Mesh1Var/Mesh1Med) },
           "dc.bin"={ c(0.9*mu.sd,1.1*mu.sd,0.5) },
           "dc.bil"={ c(log(Mesh1Med)-0.1,CV,log(Mesh1Med)+0.1,CV,0.5) },
           stop("Please provide a value of x0 (initial parameter values")
    ) }
  }


#===============================================================================
#' Calculate constant term in log-likelihood (including for non-integer counts)
#' @description Returns the constant term for binomial or multinomial counts.
#' The default is to return the logged constant.
#' The counts do not need to be integer valued.
#' @param O Matrix of observed counts
#' @export
mchoose=function(O,log=TRUE){
  O=as.matrix(O)
  if(ncol(O)==1) stop("Data must have at least 2 columns")
  rowN=apply(O,1,sum,na.rm=T)
  const=0
  for(i in 1:nrow(O)) const=const+lgamma(rowN[i]+1)-sum(lgamma(O[i,]+1),na.rm=T)
  if (log) return(const)
  else return(exp(const))
}


#===============================================================================
#' Extended binomial density function that also works for non-integer counts
#' @export
dBinom=function(y,n,p,log=FALSE) {
  if(any(y<0) | any(y>n)) stop("Infeasible data, y<0 or y>n")
  log.const=lgamma(n+1)-lgamma(y+1)-lgamma(n-y+1)
  log.prob=log.const+y*log(p)+(n-y)*log(1-p)
  log.prob[(y==0&p==0)|(y==n&p==1)]=0
  if (log) return(log.prob)
  else return(exp(log.prob))
}


#===============================================================================
#' delta method function
#' @description Returns variance matrix of functions of estimated values
#' e.g., for std errors of L50 and SR from a and b parameterization
deltamethod=function (g, mean, cov, ses = TRUE)
{
  cov <- as.matrix(cov)
  n <- length(mean)
  if (!is.list(g))
    g <- list(g)
  if ((dim(cov)[1] != n) || (dim(cov)[2] != n))
    stop(paste("Covariances should be a ", n, " by ",
               n, " matrix"))
  syms <- paste("x", 1:n, sep = "")
  for (i in 1:n) assign(syms[i], mean[i])
  gdashmu <- t(sapply(g, function(form) {
    as.numeric(attr(eval(deriv(form, syms)), "gradient"))
  }))
  new.covar <- gdashmu %*% cov %*% t(gdashmu)
  if (ses) {
    new.se <- sqrt(diag(new.covar))
    new.se
  }
  else new.covar
}


#===============================================================================
#' Modification of deltamethod (msm) that also returns function value
#' @export
#'
delta.method=function (g, mean, cov, ses = TRUE)
{
  cov <- as.matrix(cov)
  n <- length(mean)
  if (!is.list(g))
    g <- list(g)
  if ((dim(cov)[1] != n) || (dim(cov)[2] != n))
    stop(paste("Covariances should be a ", n, " by ",
               n, " matrix"))
  syms <- paste("x", 1:n, sep = "")
  for (i in 1:n) assign(syms[i], mean[i])
  gval = t(sapply(g, function(form) { as.numeric(eval(deriv(form, syms))) } ))
  gdashmu <- t(sapply(g, function(form) {
    as.numeric(attr(eval(deriv(form, syms)), "gradient")) }))
  new.covar <- gdashmu %*% cov %*% t(gdashmu)
  if (ses) { list(value=as.vector(gval),se=sqrt(diag(new.covar))) }
  else list(value=as.vector(gval),covar=new.covar)
}


#===============================================================================
#' Convert raw counts by haul to totals, and adjusts for any subsampling
#' @description Adjusts for any subsampling and (by default) sum over tows.
#' @param data Catch data in SELECT format
#' @param var.names Character vector giving the variable names.
#' The first name is always the length variable,
#' followed by names of the raw catch variables. E.g., c("lgth","nA","nB")
#' @param q.name Character vector giving the names of sampling fractions, if any.
#' @param sumHauls Sum the (adjusted) catch over all hauls?
#' @param q.ODadjust Over-dispersion adjustment. If TRUE (default) and the catches
#' have been adjusted upwards due to sampling, then the summed catch numbers are
#' adjusted downwards to have total equal to the actual number of fish measured.
#' This does not change catch share proportions, but will reduce potential
#' overfitting of polynomial and spline curves due to overdispersion.
#' @return dataframe
#' @export
Raw2Tots=function(data,var.names,q.names=NULL,sumHauls=TRUE,q.ODadjust=TRUE) {
  #For unpaired data q may be zero if count is 0. Need to avoid divide by zero
  if(!is.null(q.names)) {
    n.raw=sum(data[,var.names[-1]])
    Q=data[,q.names]
    Q[Q==0]=1
    data[,var.names[-1]]=data[,var.names[-1]]/Q
    if(q.ODadjust) {
      n.tot=sum(data[,var.names[-1]])
      data[,var.names[-1]]=data[,var.names[-1]]*n.raw/n.tot }
    }
  if(sumHauls)
    data=data[,var.names] %>% group_by(across(all_of(var.names[1]))) %>%
      summarize(across(all_of(var.names[-1]),sum)) %>% data.frame()
  data
}

#===============================================================================
#' Change long format to SELECT format
#' @description Change long format to SELECT format
#' @param by Character vector with names of the variables to join by,
#' typically TowID and lgth
#' @param gear Character giving gear variable name
#' @param freq Character giving frequency variable name
#' @param q.name Name of sampling fraction variable
#' @param paired Are data paired within gear?
#' @return dataframe, in SELECT format
#' @export
SELECT.FORMAT=function(Df,by=c("haul","lgth"),gear="gear",freq="freq",q.name=NULL,
                       paired=T) {
  UniqueCheck=Df|>group_by(across(all_of(c(by,gear))))|> summarize(m=n(),.groups="keep")
  if(max(UniqueCheck$m)>1) {
     warning("Possible data error: There are multiple rows for a unique \n combination of ",
      c(by,gear),". A unique row ID variable has been added \n")
     Df$uniqueRowID=1:nrow(Df)
     by=c("uniqueRowID",by) }
  if(is.null(q.name)) values=freq else values=c(freq,q.name)
  Wk=Df |> select(all_of( c(by,gear,values) ))
  Wk = Wk |>   pivot_wider(names_from=all_of(gear), #names_prefix=namePrefix,
                      values_from=all_of(values), values_fill=0, names_sep="")
  if(!paired) Wk[,gear]=Df[,gear]
  data.frame(Wk)
}


