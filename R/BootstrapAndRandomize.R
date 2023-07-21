## Bootstrap and randomization functions

#' Implement bootstrap
#' @description Apply the double bootstrap and return a statistic from each
#' bootstrap sample in a matrix.
#' @export
boot.SELECT=function(Data,statistic,N,SetID="Haul",BlockID=NULL,GearID=NULL,
                     Freqs=c("nfine","nwide"),within.resamp=T,...) {
  z=try( statistic(Data,...) )
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  #cat("Raw data output:",z,"\n")
  BootMatrix=matrix(NA,nrow=N,ncol=length(z))
  cat(paste("Bootstrap data: Starting a",N,"resamples bootstrap...\n"))
  if(is.null(SetID)) stop("SetID (set id variable) is required.")
  PBar <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {
    if(i%%5==0) setTxtProgressBar(PBar, i)
    bootData=Dble.boot(Data,SetID,BlockID,GearID,Freqs,within.resamp)$bootData
    boot.stat=try( statistic(bootData,...) )
    if(class(boot.stat)[1]!="try-error") BootMatrix[i,]=boot.stat
  }
  close(PBar)
  cat("\nBootstrap successfully completed\n")
  if(any(is.na(BootMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(BootMatrix)
}



#' Crude weighted average (0.25,0.5,0.25) with immediate neighbours
WgtAvg=function(y,w=c(0.25,0.5,0.25)) {
  n=length(y)
  y.left=c(y[1],y)
  y.right=c(y,y[n])
  wgt.y=w[2]*y+w[1]*y.left[1:n]+w[3]*y.right[2:(n+1)]
  wgt.y
}


#' Return a dataframe containing double bootstrapped raw data.
#' @description Double bootstrap function used by boot.SELECT function
#'
#' @param Data Stacked matrix or dataframe of catches in SELECT format,
#' with lengthclass in first column. Remaining columns are raw catch frequencies
#' @param SetID Name of grouping variable containing the tow/haul/set id number.
#' Vector of same length as number of rows in Data.
#' @param BlockID If specified, name of blocking variable.
#' Bootstrapping is first done over blocks, and then within each block.
#' @param GearID If specified, name of the gear indicator variable.
#' This restricts bootstrapping to within each gear and is intended to be used
#' for non-paired data.
#' @param within.resamp If F, no resampling at observation level.
#' @param smooth Smooth at within-haul phase to avoid losing degrees of freedom (from increasing number of zero freqs)
#'
#' @return Dataframe of double-bootstrapped freqs
#' @export
#'

Dble.boot=function(Data,SetID="Haul",BlockID=NULL,GearID=NULL,
                   Freqs=c("nfine","nwide"),within.resamp=T,smooth=F) {
  if(smooth) w=c(1,2,1)/4
  Sets=Data[,SetID]
  uniqSets=unique(Sets)
  nSets=length(uniqSets)
  #cat("\n",nSets,"sets to be double bootstrapped")
  RanList=list(uniqSets)

  if(is.null(BlockID)&is.null(GearID)) BootID=sample(uniqSets,nSets,replace=T)
  if(!is.null(BlockID)) { #Bootstraps over BlockID
    BlockSets=tapply(Sets,Data[,BlockID],unique)
    nBlocks=length(BlockSets)
    #RanList=lapply(BlockSets,sample,replace=T)
    RanList=BlockSets[sample(1:nBlocks,replace=T)]
    if(is.null(GearID)) RanList=lapply(RanList,sample,replace=T)
    BootID=unlist(RanList) }
  if(!is.null(GearID)) { #Bootstrap within Gear
    GearSets=tapply(Sets,Data[,GearID],unique)
    nGears=length(GearSets)
    WkList=as.list(1:nGears)
    for(j in 1:nGears) {
      WkList[[j]]=lapply(RanList,function(x) x[x %in% GearSets[[j]] ])
      WkList[[j]]=lapply(WkList[[j]],sample,replace=T)
    }
    BootID=unlist(WkList)
  }

  nRanSets=length(BootID)
  BootList=as.list(1:nRanSets)

  for(j in 1:nRanSets) BootList[[j]] <- Data %>% filter(Sets==BootID[j])
  if(within.resamp&smooth) {
    #Smooth to reduce number of zero obs
    #<Inefficient code since WgtAvg only needs to be done once>
    for(j in 1:nRanSets) {
      m=nrow(BootList[[j]])
      for(k in 1:length(Freqs))
        BootList[[j]][,Freqs[k]]=rpois(m,WgtAvg(BootList[[j]][,Freqs[k]],w=w))
    }
  }
  bootData=as.data.frame(rbindlist(BootList))
  if(within.resamp&!smooth) {
    m=nrow(bootData)
    for(k in 1:length(Freqs)) bootData[,Freqs[k]]=rpois(m,bootData[,Freqs[k]])
  }
  return(list(bootData=bootData,BootID=BootID))
}

#' Produce the bootstrap plot
#' @description ggplot showing fitted curve and bootstrap bounds
#'
#' @param BootPreds Matrix with bootstrap by row and fitted values at length in columns.
#' @param lenseqs Lengths at which fitted values are calculated.
#' @param predn Fitted curve
#'
#' @return ggplot GROB
#' @export
#'

BootPlot=function(BootPreds,lenseq,predn,Data=NULL,eps=0.025) {
  txt=6
  Preds.lower=apply(BootPreds,2,quantile,prob=eps,na.rm=T)
  Preds.upper=apply(BootPreds,2,quantile,prob=1-eps,na.rm=T)
  Pdf=data.frame(len=lenseq,pred=predn,low=Preds.lower,upp=Preds.upper)
  #NB geom_ribbon requires the df to have a y variable, although not used.
  BootGROB=ggplot(data=Pdf,aes(len))+
    #geom_point(data=subset(Tots,subset=c(!is.na(y))),aes(x=lgth,y=y))+
    geom_line(data=Pdf,aes(len,pred))+ylim(0,1)+
    geom_ribbon(data=Pdf,aes(x=len,ymin=low,ymax=upp),alpha=0.2)+
    xlab("Length (cm)")+ylab("Retention probability")+theme_bw()+
    theme(axis.text=element_text(size=txt),axis.title=element_text(size=txt))+
    theme(plot.margin = unit(c(0.75,0.5,0.25,0.5), "cm"))
  if(!is.null(Data)) BootGROB = BootGROB + geom_point(data=Data,aes(x=lgth,y=y))
  BootGROB
}



#' Return a permuted SELECT format dataframe
#' @description Randomization of gear type within each haul
#'
#' @param catch Stacked matrix or dataframe of catches in SELECT format
#' @param varnames Length-3 vector containing name of haulID variable, and 2 catch freq vars
#'
#' @return Dataframe, with randomized gear treatment within each haul
#' @export
#'

Randomize=function(catch,varnames=c("Haul","nwide","nfine")) {
  Tow=catch[,varnames[1]]
  uniqTows=unique(Tow)
  nTows=length(uniqTows)
  #cat("\n",nTows,"hauls to be randomized")
  RanCatchList=as.list(1:nTows)
  for(j in 1:nTows) {
    RanCatchList[[j]] = catch %>% filter(Tow==uniqTows[j])
    RanCatchList[[j]][,varnames[2:3]] <- RanCatchList[[j]][,varnames[sample(2:3)]]
  }
  ran.catch=as.data.frame(rbindlist(RanCatchList))
  return(as.data.frame(ran.catch))
}


#' Return a dataframe with permutations of haul treatment
#' @description Randomization of haul treatments (may be more than 2)
#'
#' @param Data Stacked matrix or dataframe of catches in long format. Data are not
#' assumed paired, so could be different length range within each haul.
#' @param SetID Name of grouping variable containing the tow/haul/set id number.
#' @param GearID Name of the gear (treatment) indicator variable.
#'
#' @return Dataframe, with randomized gear treatment.
#' @export
#'

RandomizeLong=function(Data,SetID="Haul",Gear="Gear") {
  Wk=Data[,c(SetID,Gear)]
  Wk[,SetID]=factor(Wk[,SetID],levels=unique(Wk[,SetID])) #Preserve order
  names(Wk)=c("SetID","Gear")
  PerSet=Wk %>% group_by(SetID,Gear) %>% summarize(n=n()) %>% data.frame()
  RanPerSet=PerSet %>% mutate(Gear=sample(Gear))
  RanGear=with(RanPerSet,rep(Gear,n))
  Data[,Gear]=RanGear
  return(Data)
}


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

#' Old workhorse
#' @export
#'


Dble.boot230425=function(Data,SetID="Haul",BlockID=NULL,
                         Freqs=c("nfine","nwide"),inner.resamp=T,smooth=F) {
  if(smooth) w=c(1,2,1)/4
  Tow=Data[,SetID]
  uniqTows=unique(Tow)
  nTows=length(uniqTows)
  #cat("\n",nTows,"hauls to be double bootstrapped")
  BootList=as.list(1:nTows)
  if(is.null(BlockID)) BootID=sample(uniqTows,nTows,replace=T)
  if(!is.null(BlockID)) {
    Block=Data[,BlockID]
    BlockTows=tapply(Tow,Block,unique)
    BootByBlock=lapply(BlockTows,sample,replace=T)
    BootID=unlist(BootByBlock)
  }
  for(j in 1:nTows) BootList[[j]] <- Data %>% filter(Tow==BootID[j])
  if(inner.resamp&smooth) {
    #Smooth to reduce number of zero obs
    #<Inefficient code since WgtAvg only needs to be done once>
    for(j in 1:nTows) {
      m=nrow(BootList[[j]])
      for(k in 1:length(Freqs))
        BootList[[j]][,Freqs[k]]=rpois(m,WgtAvg(BootList[[j]][,Freqs[k]],w=w))
    }
  }
  bootData=as.data.frame(rbindlist(BootList))
  if(inner.resamp&!smooth) {
    m=nrow(bootData)
    for(k in 1:length(Freqs)) bootData[,Freqs[k]]=rpois(m,bootData[,Freqs[k]])
  }
  return(list(bootData=bootData,BootID=BootID))
}
