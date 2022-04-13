# Paused at development of boot.SELGAM. Search Fawk

#' Fit selection curves to fishing gears using the SELECT method.
#' @description SELECT is an acronym for the Share Each LEngths Catch Total method
#' of Millar (1992). SELECT is a general statistical framework to  estimate the
#' selectivity of trawls, hooks and gillnets from experimental data.
#' @references Millar, R. B. (1992).
#' Estimating the size-selectivity of fishing gear by conditioning on the
#' total catch. Journal of the American Statistical Association. 87: 962-968.
#' @export
SELECT=function(Data,dtype,stype="logistic",Q=NULL,Meshsize=NULL,x0=NULL,rel.power=NULL,
                penalty.func=NULL,print.out=T,
                control=list(maxit=1000,reltol=1e-8),fit=T) {
  #Assemble the SELECT arguments, and put in SELECT.args
  dtype=substr(dtype,1,2)
  rtype=paste0(dtype,".",stype)
  r=propncurves(rtype) #Get propn catch curve function
  if(!is.null(Q)) Data[,-1]=Data[,-1]/Q
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
  if(ncol(Counts)!=length(Meshsize))
    stop("Number of mesh sizes should be ",ncol(Counts))
  if(is.null(x0)) x0=StartValues(rtype,Data)
  SELECT.args=list(Data=Data,rtype=rtype,Q=Q,Meshsize=Meshsize,x0=x0,
                   rel.power=rel.power,penalty.func=penalty.func)
  ll.init=Const-nllhood(theta=x0,Data,Meshsize,r,rel.power,penalty.func)
  #If not fitting (using just x0)
  if(!fit) {
    if(print.out) cat("SELECT model evaluated at x0 - no optimization: \n")
    z=c(list(par=x0),SELECT.args,loglik.init=ll.init)
    return( invisible(z) ) }
  #Calculate logliks of full model and at x0
  CountTotals=apply(Counts,1,sum,na.rm=TRUE)
  CountTotals=ifelse(CountTotals==0,Inf,CountTotals)
  CountPropns=Counts/CountTotals
  fullfit.l=Const+sum(Counts[Counts>0]*log(CountPropns[Counts>0]),na.rm=TRUE)
  if(print.out) cat("Evaluated at initial parameters",round(x0,2),
                    ", the log-likelihood is",ll.init,"\n")
  #Fit=TRUE
  if(print.out)
    cat("Fitting SELECT model with",stype,"selection curves to",design,"data.\n")
  fit=optim(x0,nllhood,Data=Data,Meshsize=Meshsize,r=r,rel.power=rel.power,
            penalty.func=penalty.func,hessian=T,control=control)
  fit$value=fit$value-Const
  fit=fit[names(fit) %in% c("counts", "message") == FALSE]
  Dev=2*(fullfit.l+fit$value)
  if(print.out) {
    cat(paste0("Convergence code ",fit$convergence,
             ": Optimizer has ",ifelse(fit$convergence==0,"","*NOT*"),"converged\n"))
    cat("Pars=",fit$par,", Deviance=",Dev,", #len classes=",sum(CountTotals>0),"\n") }
  z=c(Call=match.call(),fit=fit,deviance=Dev,SELECT.args,loglik.init=ll.init)
  class(z)="SELECT"
  return(invisible(z))
#  invisible(c(fit,deviance=Dev,rtype=rtype,rel.power=list(rel.power),
#              Meshsize=list(Meshsize),Data=list(Data)))
#  invisible(list(fit=fit,deviance=Dev,Call=SELECTcall,Data=Data,rtype=rtype,
#              rel.power=rel.power,Meshsize=Meshsize) )
}
#===============================================================================

#' @export
boot.SELECT=function(Data,statistic,N,SetID="Haul",Freqs=c("nfine","nwide")) {
  z=try( statistic(Data) )
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  #cat("Raw data output:",z,"\n")
  BootMatrix=matrix(NA,nrow=N,ncol=length(z))
  cat(paste("Bootstrap data: Starting a",N,"resamples bootstrap...\n"))
  if(is.null(SetID)) stop("SetID (set id variable) is required.")
  PBar <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {
    if(i%%5==0) setTxtProgressBar(PBar, i)
    bootData=Dble.boot(Data,SetID,Freqs,smooth=F)$bootData
    boot.stat=try( statistic(bootData) )
    if(class(boot.stat)[1]!="try-error") BootMatrix[i,]=boot.stat
  }
  close(PBar)
  cat("\nBootstrap successfully completed\n")
  if(any(is.na(BootMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(BootMatrix)
}



boot.SELECTold=function(Data,dtype,stype="logistic",SetID=NULL,nboots=1000,
                     Q=NULL,Meshsize=NULL,
                     x0=NULL,rel.power=NULL,penalty.func=NULL,print.out=F) {
  cat("Raw data:")
  z=SELECT(Data,dtype,stype,Q=Q,Meshsize=Meshsize,
           x0=x0,rel.power=rel.power,penalty.func=penalty.func,print.out=T)
  cat(paste("Bootstrap data: Starting a",nboots,"resamples bootstrap...\n"))
  if(is.null(SetID)) stop("SetID (set id variable) is required.")
  Est=Estimates(z)
  nEst=nrow(Est)
  BootEstMatrix=matrix(NA,nrow=nboots,ncol=nEst)
  colnames(BootEstMatrix)=rownames(Est)
  PBar <- txtProgressBar(min = 0, max = nboots, style = 3)
  for(i in 1:nboots) {
    if(i%%5==0) setTxtProgressBar(PBar, i)
    Boot.df=Dble.boot(Data,SetID,smooth=F)$bootData
    if(!is.null(Q)) Boot.df[,-1]=Boot.df[,-1]/Q #Fawk. Q not booted
    #Aggregate for speed
    colnames(Boot.df)[1]="size"
    Agg.tib <- Boot.df %>% group_by(size) %>% summarize_all(sum)
    Agg.df=data.frame(Agg.tib)
    boot.Est=try( Estimates( SELECT(Agg.df,dtype,stype,Meshsize=Meshsize,x0=x0,
                rel.power=rel.power, penalty.func=penalty.func,print=F) ),
                silent=T )
    if(class(boot.Est)[1]!="try-error") BootEstMatrix[i,]=boot.Est[,1]
  }
  close(PBar)
  cat("\nBootstrap successfully completed\n")
  if(any(is.na(BootEstMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(BootEstMatrix)
}
  #if(class(boot.fit)[1]!="try-error") BootEstMatrix[i,]=Estimates(boot.fit)[,1]


boot.SELGAM=function(nboots=1000,SetID=NULL,len.seq=NULL,data,Q=NULL,...) {
  Data=data
  if(!is.null(Q)) Data[,1]=Data[,1]/Q
  z=try( gam(data=Data,...) ) #Not yet corrected for Q
  if(class(z)[1]=="try-error") stop("GAM fit to actual data was unsuccessful")
  cat(paste("Bootstrap data: Starting a",nboots,"resamples bootstrap...\n"))
  if(is.null(SetID)) stop("SetID (set id variable) is required.")
  if(is.null(len.seq)) stop("Length vector for prediction is required")
  BootPredMatrix=matrix(NA,nrow=nboots,ncol=length(len.seq))
  #colnames(BootEstMatrix)=rownames(Est)
  PBar <- txtProgressBar(min = 0, max = nboots, style = 3)
  for(i in 1:nboots) {
    if(i%%5==0) setTxtProgressBar(PBar, i)
    Boot.df=Dble.boot(data,SetID,smooth=F)$bootData
    if(!is.null(Q)) Boot.df[,-1]=Boot.df[,-1]/Q
    #Aggregate for speed
    #colnames(Boot.df)[1]="size"
    Agg.tib <- Boot.df %>% group_by(size) %>% summarize_all(sum)
    Agg=data.frame(Agg)
    Agg$n=Agg[,2]+Agg[,3]; Agg$y=Agg[,3]/Agg$n;
    boot.Est=try( predict( gam(data=Agg.df,...),type="response"), silent=F )
    if(class(boot.Est)[1]!="try-error") BootEstMatrix[i,]=boot.Est
  }
  close(PBar)
  cat("\nBootstrap successfully completed\n")
  if(any(is.na(BootPredMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(BootPredMatrix)
}
#boot.SELGAM(y~s(length,bs="bs",m=c(2,1)),family=quasibinomial,weights=n,data=Totboot.df)

WgtAvg=function(y,w=c(0.25,0.5,0.25)) {
  n=length(y)
  y.left=c(y[1],y)
  y.right=c(y,y[n])
  wgt.y=w[2]*y+w[1]*y.left[1:n]+w[3]*y.right[2:(n+1)]
  wgt.y
}
#===============================================================================



#print.SELECT=function(z,...) {
#  keep=c("Call","rtype","par","convergence","hessian","deviance")
#  zsub=z[names(z) %in% keep]
#  print.default(zsub)
#}


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
  nll=-sum(Counts[Counts>0]*log(phi[Counts>0]))
  #nll=sum(byrow_dmultinom(Counts,prob=phi,log=T))
  nll=nll+penalty.func(theta)
  return(nll) }


#' Provide a diagnostic summary of the SELECT model fit
#' @description Provides log-likelihoods, deviances,
#' over-dispersion correction factors (if minE>0), and a plot of deviance residuals.
#' @export
ModelCheck=function(fit,minE=0,xlab="Length (cm)",ylab = "Propn in exptl gear",
                    print.out=T,plots=T,pex=1,...) {
  r=propncurves(fit$rtype) #Get propn catch curve function
  lens=fit$Data[,1]; nlens=length(lens)
  Meshsize=fit$Meshsize; nmeshes=length(Meshsize)
  O=fit$Data[,-1]; #Matrix of observed counts
  Const=mchoose(O,log=T) #Constant for log-likelihoods
  Ototals=apply(O,1,sum,na.rm=TRUE)
  Opropns=O/Ototals
  NullPropns=matrix(apply(O,2,sum,na.rm=TRUE),nrow=nlens,ncol=nmeshes,byrow=T)/sum(O)
  rmatrix=outer(lens,Meshsize,r,fit$par)
  rmatrix[is.na(O)]=NA #No fitted retention for missing meshsizes
  rmatrix=t(t(rmatrix)*fit$rel.power)
  phi=rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
  E=apply(O,1,sum,na.rm=TRUE)*phi #Matrix of expected counts
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
  dof=nrow(NonZeroDat)*(nmeshes-1)-length(fit$par)-sum(is.na(NonZeroDat))
  Deviance.CF=Deviance/dof; Pearson.CF=Pearson.chisq/dof
  out1=cbind(null.l,model.l,full.l,Deviance,Pearson.chisq,dof,Deviance.CF,Pearson.CF)
  outlist=list(stats=out1,fit=E)
  #If n cells for a given length have freq>=minE, it contributes max(0,n-1) dof
  if(minE>0) {
    Index=(E>minE)
    RowIndex=(1:nrow(E))[apply(Index,1,sum,na.rm=TRUE)==1]
    Index[RowIndex,]=FALSE
    dof=sum(pmax(0,apply(Index,1,sum,na.rm=TRUE)-1))-length(fit$par)
    Pearson.chisq=sum((Pearson.resids^2)[Index],na.rm=TRUE)
    Deviance=sum((Dev.resids^2)[Index],na.rm=TRUE)
    Deviance.CF=Deviance/dof
    Pearson.CF=Pearson.chisq/dof
    out2=cbind(Deviance,Pearson.chisq,dof,Deviance.CF,Pearson.CF)
    outlist=list(stats=out1,fit=E)
  }
  #if(plots) {
  AreLensUnique=(length(lens)==length(unique(lens)))
  xyticks=c(nlens-1,10,7)
  if(nmeshes>2&AreLensUnique) {
    #Plot deviance residuals
    plot(1,1,xlim=range(lens),xlab=xlab,ylab="Deviance residuals",
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
      abline(h=0)
      #Plot the fit to propns
      Data=fit$Data; pwr=fit$rel.power
      plotlens=seq(min(lens),max(lens),length=100)
      plot(lens,Data[,3]/(Data[,2]+Data[,3]),type=ifelse(AreLensUnique,"b","p"),
           ylim=c(0,1),xlab=xlab,ylab=ylab,lab=xyticks,...)
      fit.rmatrix=outer(plotlens,fit$Meshsize,r,fit$par)
      fit.rmatrix=t(t(fit.rmatrix)*pwr)
      phi=fit.rmatrix[,2]/(fit.rmatrix[,1]+fit.rmatrix[,2])
      lines(plotlens,phi,type="l",lty=2) }
  #} #End of plots
  if(print.out) {
    cat("Model fit:\n"); print(out1[1,]);
    if(minE>0) {
      cat("\nCorrection factors from cells with expected count >",minE,":\n");
      outlist=list(stats=out1,CF=out2,fit=E)
      print(out2[1,]) } }
  invisible(outlist)
}



#' Returns selectivity parameters
#' @description Calculates the estimates and standard error of selectivity
#' curve implemented by SELECT. The user can edit to add new calculations for
#' any user-supplied retention curves.
#' @export
Estimates=function(fit) {
  x=fit$par;
  if(is.null(fit$hess))
    stop("\n SELECT fit does not include a hessian. Perhaps you used fit=F.\n")
  varx=solve(fit$hess)
  names=c("Mode(mesh1)","Std dev.(mesh1)")
  rtype=fit$rtype
  switch( #"re" and "cc" estimates have identical parameter sets
    {rtype=ifelse(rtype=="cc.logistic","re.logistic",rtype)
     rtype=ifelse(rtype=="cc.richards","re.richards",rtype)},
    "re.norm.loc"={ pars=x; varpars=varx },
    "re.norm.sca"={ pars=x; varpars=varx },
    "re.gamma"={
      pars=c((x[1]-1)*x[2],sqrt(x[1]*x[2]^2))
      varpars=deltamethod(list(~(x1-1)*x2,~x2*sqrt(x1)),x,varx,ses=F)},
    "re.logistic"={ pars=x; varpars=varx; names=c("L50","SR") },
    "re.lognorm"={
      pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)))
      varpars=deltamethod(list(~exp(x1-x2^2),
               ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),x,varx,ses=F)},
    "re.binorm.sca"={
      pars=c(x[1:4],exp(x[5])/(1+exp(x[5])))
      names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                    "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
      varpars=deltamethod(list(~x1,~x2,~x3,~x4,~exp(x5)/(1+exp(x5))),
                                                x,varx,ses=F)},
    "re.bilognorm"={
      pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)),
             exp(x[3]-x[4]^2),sqrt(exp(2*x[3]+x[4]^2)*(exp(x[4]^2)-1)),
             exp(x[5])/(1+exp(x[5])))
      names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                    "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
      varpars=deltamethod(
        list(~exp(x1-x2^2),~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1)),
             ~exp(x3-x4^2),~sqrt(exp(2*x3+x4^2)*(exp(x4^2)-1)),
             ~exp(x5)/(1+exp(x5))),x,varx,ses=F)},
    "re.richards"={
      pars=c(x[1],x[2],exp(x[3]))
      names=c("L50","SR","delta")
      varpars=deltamethod(list(~x1,~x2,~exp(x3)),x,varx,ses=F)},
    "ph.logistic"={
      pars=c(x[1],x[2],plogis(x[3]))
      names=c("L50","SR","p")
      varpars=deltamethod(list(~x1,~x2,~exp(x3)/(1+exp(x3))),x,varx,ses=F)},
    "ph.richards"={
      pars=c(x[1],x[2],exp(x[3]),plogis(x[4]))
      print(pars)
      names=c("L50","SR","delta","p")
      varpars=deltamethod(list(~x1,~x2,~exp(x3),~exp(x4)/(1+exp(x4))),x,varx,ses=F)},
    ###################Deprecated options#######################################
    "re.logistic.ab"={ #If using a,b parameterization (old code)
      pars=c(-x[1]/x[2],2*(log(3))/x[2])
      names=c("L50","SR")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2),x,varx,ses=F)},
    "re.richards.ab"={ #If using a,b parameterization (old code)
      delta=exp(x[3])
      pars=c((log(0.5^delta/(1-0.5^delta))-x[1])/x[2],
             (log(0.75^delta/(1-0.75^delta))-log(0.25^delta/(1-0.25^delta)))/x[2],
             delta)
      names=c("L50","SR","delta")
      varpars=deltamethod(list(
        ~(log(0.5^exp(x3)/(1-0.5^exp(x3)))-x1)/x2,
        ~(log(0.75^exp(x3)/(1-0.75^exp(x3))))/x2
        -(log(0.25^exp(x3)/(1-0.25^exp(x3))))/x2,~exp(x3)),x,varx,ses=F)},
    "ph.logistic.ab"={
      pars=c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
      names=c("L50","SR","p")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
                          x,varx,ses=F)},
    "ph.richards.ab"={
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
      'Possible relative selection types are , "norm.loc", "norm.sca", "gamma",
        "lognorm", "logistic", "richards, "binorm.sca", and "bilognorm". \n',
      'Possbile covered-codend and alternative hauls types are
                   "logistic" and "richards" \n'))
  )#End of switch
  estimates=cbind(pars,sqrt(diag(varpars)))
  colnames(estimates)=c("par","s.e.")
  rownames(estimates)=names
  return(estimates) }



#' Plot the fitted selectivity curve(s)
#' @export
PlotCurves=function(fit,plotlens=NULL,Meshsize=NULL,rel.power=NULL,standardize=F,
                    plot.out=T,xlab="Length (cm)",ylab="Retention curve",...) {
  s=selncurves(fit$rtype) #Get selection curve function
  if(is.null(plotlens)) plotlens=fit$Data[,1]
  if(is.null(Meshsize)) Meshsize=fit$Meshsize
  if(is.null(rel.power)) pwr=fit$rel.power
  smatrix=outer(plotlens,Meshsize,s,fit$par)
  smatrix=t(t(smatrix)*pwr)
  if(standardize) smatrix=smatrix/max(smatrix)
  #Plot propn retained if only two gears
  if(plot.out){
    if(length(Meshsize)==2) {
      plot(plotlens,smatrix[,2],type="l",ylim=c(0,1),xlab=xlab,ylab=ylab,...)
      abline(h=c(0.25,0.5,0.75),lty=3) }
    else {
      matplot(plotlens,smatrix,type="l",ylim=c(0,1),xlab=xlab,ylab=ylab,...) }
  }
  lensmatrix=cbind(plotlens,smatrix)
  colnames(lensmatrix)=c("Length",paste0("Gear",1:length(Meshsize)))
  invisible(lensmatrix) }


#' @export
StartValues=function(rtype,Data) {
  switch(substr(rtype,1,6),
         #"logist"={ c(mean(rep(Data[,1],Data[,3])),sd(rep(Data[,1],Data[,3]))) },
         "cc.log"={ c(mean(rep(Data[,1],Data[,3])),sd(rep(Data[,1],Data[,3]))) },
         "cc.ric"={ c(mean(rep(Data[,1],Data[,3])),sd(rep(Data[,1],Data[,3])),0) },
         "ph.log"={ c(mean(rep(Data[,1],Data[,3])),sd(rep(Data[,1],Data[,3])),0) },
         "ph.ric"={ c(mean(rep(Data[,1],Data[,3])),sd(rep(Data[,1],Data[,3])),0,0) },
         stop("Please provide a value of x0 (initial parameter values")
  )
}

#' @export
mchoose=function(x,log=TRUE){
  x=as.matrix(x)
  if(ncol(x)==1) stop("x[] must have at least 2 columns")
  rowN=apply(x,1,sum)
  const=0
  for(i in 1:nrow(x)) const=const+lgamma(rowN[i]+1)-sum(lgamma(x[i,]+1))
  if (log) return(const)
  else return(exp(const))
}

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

#' @export
byrow_dmultinom=function(x, prob, log = FALSE) {
  x=as.matrix(x); prob=as.matrix(prob)
  if(nrow(x)!=nrow(prob)|ncol(x)!=ncol(prob))
    stop("x[] and prob[] must have same dimension.")
  if(ncol(x)==1) stop("x[] must have at least 2 columns")
    k=ncol(x)
    rowN=apply(x,1,sum)
    r=sapply(1:nrow(x), function(i) {
      xx <- x[i,]
      pp <- prob[i, ]
      i0 <- pp == 0
      if (any(i0)) {
        if (any(xx[i0] != 0))
          return(0)
        xx <- xx[!i0]
        pp <- pp[!i0]
      }
      return(lgamma(rowN[i] + 1) + sum(xx*log(pp) - lgamma(xx+1)))
    })
    if (log)
      return(r)
    else return(exp(r))
}
#fullfit.l=sum(byrow_dmultinom(Counts,prob=CountPropns,log=TRUE))


#' @export
dmultinomial=function (x, size = NULL, prob, log = FALSE)
{
  if (length(x) == 0)
    return(numeric(0))
  if (is.vector(prob))
    prob <- t(prob)
  if (is.vector(x))
    x <- t(x)
  K <- ncol(prob)
  maxn <- max(nx <- nrow(x), np <- nrow(prob))
  if (ncol(x) != K)
    stop("x[] and prob[] must be equal length vectors or equal col matrix.")
  if (nx != maxn)
    x <- matrix(t(x), ncol = K, nrow = maxn, byrow = TRUE)
  if (np != maxn)
    prob <- matrix(t(prob), ncol = K, nrow = maxn, byrow = TRUE)
  N <- apply(x, 1, sum)
  if (is.null(size))
    size <- N
  size <- rep(size, length = maxn)
  if (any(size != N))
    stop("at least one size != sum(x), i.e. one is wrong")
  if (any(prob < 0) || any((s <- apply(prob, 1, sum)) == 0))
    stop("probabilities cannot be negative nor all 0.")
  prob <- prob/s
  if (any(as.integer(x + 0.5) < 0))
    stop("'x' must be non-negative")
  r <- sapply(1:maxn, function(y) {
    xx <- as.integer(x[y, ] + 0.5)
    pp <- prob[y, ]
    i0 <- pp == 0
    if (any(i0)) {
      if (any(xx[i0] != 0))
        return(0)
      xx <- xx[!i0]
      pp <- pp[!i0]
    }
    return(lgamma(size[y] + 1) + sum(xx * log(pp) - lgamma(xx + 1)))
  })
  if (log)
    return(r)
  else return(exp(r))
}

