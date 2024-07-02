## Model evaluation and helper functions
# ModelCheck, Estimates, PlotCurves
# nllhood, calcOD, StartValues, mchoose, deltamethod, Raw2Tots & SELECT_FORMAT

#===============================================================================
#' Provide a diagnostic summary of the SELECT model fit
#' @description Provides log-likelihoods, deviances,
#' over-dispersion correction factors (if minE>0), and a plot of deviance residuals.
#' @export
ModelCheck=function(fit,minE=0,xlab="Length (cm)",ylab = "Propn in exptl gear",
                    print.out=T,plots=c(T,T),plotlens=NULL,pex=1,...) {
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
    #But also need to exclude a cell if it is the only one in a row with freq>=minE
    RowIndex=(1:nrow(E))[apply(Index,1,sum,na.rm=TRUE)==1]
    Index[RowIndex,]=FALSE
    dof=sum(pmax(0,apply(Index,1,sum,na.rm=TRUE)-1))-length(fit$par)
    Pearson.chisq=sum((Pearson.resids^2)[Index],na.rm=TRUE)
    Deviance=sum((Dev.resids^2)[Index],na.rm=TRUE)
    Deviance.CF=Deviance/dof
    Pearson.CF=Pearson.chisq/dof
    out2=cbind(Deviance,Pearson.chisq,dof,Deviance.CF,Pearson.CF)
    #outlist=list(stats=out1,fit=E)
  }
  #plots argument controls plotting of both deviance and fits

  if(length(plots)==1) plots=c(plots,plots)
  xyticks=c(nlens-1,10,7)
  AreLensUnique=(length(lens)==length(unique(lens)))
  if(plots[1]) { #Plot deviance residuals
  if(nmeshes>2&AreLensUnique) {
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
    cat("Model fit:\n"); print(out1[1,]);
    if(minE>0) {
      cat("\nCorrection factors from cells with expected count >",minE,":\n")
      print(out2[1,])
	    outlist=c(outlist,CF=out2) } }
  invisible(outlist)
}


#===============================================================================
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
	"ph.logistic"={
      pars=c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
      names=c("L50","SR","p")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
                          x,varx,ses=F)},
	"re.logistic"={
      pars=c(-x[1]/x[2],2*(log(3))/x[2])
      names=c("L50","SR")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2),x,varx,ses=F)},
  "ph.richards"={
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
	"re.richards"={
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
	"ph.logistic.L50SR"={
      pars=c(x[1],x[2],plogis(x[3]))
      names=c("L50","SR","p")
      varpars=deltamethod(list(~x1,~x2,~exp(x3)/(1+exp(x3))),x,varx,ses=F)},
    "ph.richards.L50SR"={
      pars=c(x[1],x[2],exp(x[3]),plogis(x[4]))
      print(pars)
      names=c("L50","SR","delta","p")
      varpars=deltamethod(list(~x1,~x2,~exp(x3),~exp(x4)/(1+exp(x4))),x,varx,ses=F)},
	"re.logistic.L50SR"={ pars=x; varpars=varx; names=c("L50","SR") },
    "re.logistic.ab"={ #If using a,b parameterization (old code)
      pars=c(-x[1]/x[2],2*(log(3))/x[2])
      names=c("L50","SR")
      varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2),x,varx,ses=F)},
	"re.richards.L50SR"={
      pars=c(x[1],x[2],exp(x[3]))
      names=c("L50","SR","delta")
      varpars=deltamethod(list(~x1,~x2,~exp(x3)),x,varx,ses=F)},
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



#===============================================================================
#' Plots the fitted curves and returns fitted values at specified lengths
#' @description Plot the fitted selectivity curve(s)
#' @export
PlotCurves=function(fit,plotlens=NULL,Meshsize=NULL,rel.power=NULL,standardize=F,
                    plot.out=T,xlab="Length (cm)",ylab="Retention curve",type="l",
                    ylim=c(0,1),...) {
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
      plot(plotlens,smatrix[,2],ylim=ylim,xlab=xlab,ylab=ylab,type=type,...)
      abline(h=c(0.25,0.5,0.75),lty=3) }
    else {
      matplot(plotlens,smatrix,ylim=ylim,xlab=xlab,ylab=ylab,type=type,...) }
  }
  lensmatrix=cbind(plotlens,smatrix)
  colnames(lensmatrix)=c("Length",paste0("Gear",1:length(Meshsize)))
  if(plot.out) invisible(lensmatrix)
  if(!plot.out) lensmatrix
}



#===============================================================================
#' Negative log-likelihood of the SELECT model (excluding constants)
#' @description The general negative log-likelihood function. Provided for
#' completeness, but not intended for user use.
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
#' Function to calculate overdispersion
#' @description Returns OD estimate from matrices of observed and expected counts
#' @export
calcOD=function(O,E,npar,minE=1) {
  if(!identical(dim(O),dim(E))) stop("O and E must be same dimension")
  if( max(abs(apply(O,1,sum)-apply(E,1,sum)))>1e-6 )
     cat("\nWARNING: Row sums of observed and expected not all equal\n")
  Pearson.resids=(O-E)/sqrt(E)
  wk=O*log(O/E); wk[is.na(wk)]=0
  Dev.resids=sign(O-E)*sqrt(2*(E-O+wk))
  #If n cells for a given length have freq>=minE, it contributes max(0,n-1) dof
  Index=(E>minE)
  #But also need to exclude a cell if it is the only one in a row with freq>=minE
  RowIndex=(1:nrow(E))[apply(Index,1,sum,na.rm=TRUE)==1]
  Index[RowIndex,]=FALSE
  dof=sum(pmax(0,apply(Index,1,sum,na.rm=TRUE)-1))-npar
  Pearson.chisq=sum((Pearson.resids^2)[Index],na.rm=TRUE)
  Deviance=sum((Dev.resids^2)[Index],na.rm=TRUE)
  Deviance.CF=Deviance/dof
  Pearson.CF=Pearson.chisq/dof
  ODstats=cbind(Deviance,Pearson.chisq,dof,Deviance.CF,Pearson.CF)
  return(ODstats)
  }

#===============================================================================
#' Initial parameter values for logistic and Richards curves
#' @description Returns crude data-driven starting values for logistic and
#' Richards curves
StartValues=function(rtype,Data) {
  CodendMean=mean(rep(Data[,1],Data[,3]))
  CodendSd=sd(rep(Data[,1],Data[,3]))
  b0=2*log(3)/CodendSd; a0=-b0*CodendMean
  switch(substr(rtype,1,6),
         "cc.log"={ c(a0,b0) },
         "cc.ric"={ c(a0,b0,0) },
         "ph.log"={ c(a0,b0,0) },
         "ph.ric"={ c(a0,b0,0,0) },
         stop("Please provide a value of x0 (initial parameter values")
  )
}


#===============================================================================
#' Constant term in log-likelihood
#' #' @export
mchoose=function(x,log=TRUE){
  x=as.matrix(x)
  if(ncol(x)==1) stop("x[] must have at least 2 columns")
  rowN=apply(x,1,sum,na.rm=T)
  const=0
  for(i in 1:nrow(x)) const=const+lgamma(rowN[i]+1)-sum(lgamma(x[i,]+1),na.rm=T)
  if (log) return(const)
  else return(exp(const))
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
#' @export
Raw2Tots=function(data,var.names,q.names=NULL,useTots=T) {
  Counts=data[,var.names[-1]]
  #For unpaired data q may be zero if count is 0. Need to avoid divide by zero
  if(!is.null(q.names)) {
    Q=data[,q.names]
    Q[Q==0]=1e-12
    data[,var.names[-1]]=data[,var.names[-1]]/Q }
  if(useTots) {
    Data=data[,var.names] %>% group_by(across(all_of(var.names[1]))) %>%
      summarize(across(all_of(var.names[-1]),sum)) %>% data.frame() }
  Data
}

#===============================================================================
#' Put length-freqs stacked by gear (and/or other variables) into SELECT format.
#' Default assumes input df has columns named "TowID", "lgth", "freq" and "gear"
#' where gear may be multilevel.
#' @description Change long format to SELECT format
#'
#' @param by Character vector with names of the variables to join by,
#' typically TowID and lgth
#' @param gear Character giving gear variable name
#' @param freq Character giving frequency variable name
#'
#' @return Dataframe, in SELECT format
#' @export
#'

SELECT_FORMAT240614=function(Df,by=c("TowID","lgth"),gear="gear",freq="freq") {
  wk=split(Df,Df[,gear])
  ngear=length(wk)
  nby=length(by)
  freq.names=paste0("n",names(wk))
  Stacked.df=full_join(wk[[1]][,c(by,freq)],wk[[2]][,c(by,freq)],by=by)
  if(ngear>2) {
    for(k in 3:ngear)
      Stacked.df=full_join(Stacked.df,wk[[k]][,c(by,freq)],by=by) }
  names(Stacked.df)[nby+(1:ngear)]=freq.names
  #Reshuffle columns so that length is in column 1
  Stacked.df=Stacked.df[,c(nby+(0:ngear),1:(nby-1))]
  Stacked.df[is.na(Stacked.df)] <- 0
  Stacked.df
}

SELECT_FORMAT=function(Df,by=c("haul","lgth"),gear="gear",freq="freq",q.name=NULL,
                       paired=T) {
  Wk=Df
  if(!paired) Wk=Wk |> mutate(uniqRowID=row_number())
  if(!is.null(q.name)) freq=c(freq,q.name)
  namePrefix=c("n","q")[1:length(freq)]
  Wk=Wk |> group_by(across(all_of(by))) |>
           pivot_wider(names_from=all_of(gear), #names_prefix=namePrefix,
                      values_from=all_of(freq), values_fill=0, names_sep="") |>
           select(-uniqRowID)
  if(!paired) Wk[,gear]=Df[,gear]
  Wk
}
#X=SELECT_FORMAT(Df,by=c("Haul","lgth"),gear="Gear",freq="n",q.name="q"); head(X)
