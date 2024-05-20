## Bootstrap and randomization functions

#' Bootstrap catch data
#' @description Applies a double bootstrap to data in SELECT format and
#' evaluates the vector valued function `statistic`. The returned value is a
#' nsim by length(statistic) matrix of bootstrap statistics.
#' @export
bootSELECT=function(obj,statistic,haul="Haul",nsim=2,
                     block=NULL,gear=NULL,within.resamp=TRUE,...) {
  data=obj$data
  Freqs=obj$var.names[-1] #lgth is not needed
  z=try( statistic(data,...) )
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  #cat("Raw data output:",z,"\n")
  BootMatrix=matrix(NA,nrow=nsim,ncol=length(z))
  cat(paste("\nBootstrap data: Starting a",nsim,"resamples bootstrap...\n"))
  if(is.null(haul)) stop("haul (set id variable) is required.")
  PBar <- txtProgressBar(min = 0, max = nsim, style = 3)
  for(i in 1:nsim) {
    if(i%%5==0) setTxtProgressBar(PBar, i)
    bootData=Dble.boot(data,haul,block,gear,Freqs,within.resamp)$bootData
    boot.stat=try( statistic(bootData,...) )
    #cat("Bootstrap data output:",boot.stat,"\n")
    if(class(boot.stat)[1]=="try-error") {
      cat("\nError running on bootstrap",i,"data\n")
      print(head(bootData)) }
    if(class(boot.stat)[1]!="try-error") BootMatrix[i,]=boot.stat
  }
  close(PBar)
  cat("\nBootstrap successfully completed\n")
  if(any(is.na(BootMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(BootMatrix)
}



#' Permute catch data
#' @description Permutes data in SELECT format and returns the value of the user
#' supplied function`statistic`in a nsim by length(statistic) matrix
#' @export
permSELECT=function (obj,statistic,haul="Haul", nsim=2,...)
{
  data=obj$data
  Freqs=obj$var.names[-1]
  z = try(statistic(data, ...))
  if (class(z)[1] == "try-error")
    stop("Error running statistic on actual data")
  PermMatrix = matrix(NA, nrow = nsim, ncol = length(z))
  Freqs=var.names[-1] #lgth is not needed
  cat(paste("\nStarting on", nsim, "permutations...\n"))
  if (is.null(haul))
    stop("haul is required.")
  PBar <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (i in 1:nsim) {
    if (i%%5 == 0)
      setTxtProgressBar(PBar, i)
    permData = Randomize(data, varnames=c(haul,Freqs))
    perm.stat = try(statistic(permData, ...))
    if (class(perm.stat)[1] == "try-error") {
      cat("\nError running on permutation", i, "data\n")
      print(head(permData))
    }
    if (class(perm.stat)[1] != "try-error")
      PermMatrix[i, ] = perm.stat
  }
  close(PBar)
  cat("\nPermutations successfully completed\n")
  if (any(is.na(PermMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(PermMatrix)
}



#' Crude weighted average (0.25,0.5,0.25) with immediate neighbours
WgtAvg=function(y,w=c(0.25,0.5,0.25)) {
  n=length(y)
  y.left=c(y[1],y)
  y.right=c(y,y[n])
  wgt.y=w[2]*y+w[1]*y.left[1:n]+w[3]*y.right[2:(n+1)]
  wgt.y
}

#' Avoid issue with sample from a single value
SafeSample=function(x,replace=T) {
  if(length(x)>1) x=sample(x,replace=replace)
  return(x)
}


#' Return a dataframe containing double bootstrapped raw data.
#' @description Double bootstrap function used by boot.SELECT function
#'
#' @param Data Stacked matrix or dataframe of catches in SELECT format,
#' with lengthclass in first column. Remaining columns are raw catch frequencies
#' @param haul Name of grouping variable containing the tow/haul/set id number.
#' Vector of same length as number of rows in Data.
#' @param BlockID If specified, name of blocking variable.
#' Bootstrapping is first done over blocks, and then sets within each block.
#' @param GearID If specified, name of the gear indicator variable.
#' This restricts bootstrapping of sets to within each gear and is intended to be used
#' for non-paired data.
#' @param within.resamp If F, no resampling at observation level.
#' @param smooth Smooth at within-haul phase to avoid losing degrees of freedom (from increasing number of zero freqs)
#'
#' @return Dataframe of double-bootstrapped freqs
#' @export
#'

Dble.boot=function(data,haul="Haul",BlockID=NULL,GearID=NULL,
                   Freqs=c("nfine","nwide"),within.resamp=T,smooth=F) {
  Data=as.data.frame(data) #So that code will work with a tibble
  if(smooth) w=c(1,2,1)/4
  Sets=Data[,haul]
  uniqSets=unique(Sets)
  nSets=length(uniqSets)
  #cat("\n",nSets,"sets to be double bootstrapped")
  RanList=list(uniqSets)

  #Resample the Sets, subject to GearID and BlockID where applicable
  if(is.null(BlockID)&is.null(GearID)) BootID=sample(uniqSets,nSets,replace=T)
  if(!is.null(BlockID)) { #Bootstraps over BlockID
    BlockSets=tapply(Sets,Data[,BlockID],unique)
    nBlocks=length(BlockSets)
    RanList=BlockSets[sample(1:nBlocks,replace=T)]
    #Sample the Sets within BlockID
    if(is.null(GearID)) RanList=lapply(RanList,SafeSample,replace=T)
    BootID=unlist(RanList) }
  if(!is.null(GearID)) { #Bootstrap within Gear
    GearSets=tapply(Sets,Data[,GearID],unique)
    GearSets=lapply(GearSets,SafeSample)
    BootID=unlist(GearSets)
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
  bootData=as.data.frame(bind_rows(BootList))
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

BootPlot=function(BootPreds,lenseq,predn,Data=NULL,eps=0.025,txt=8,
                  xlab="Length (cm)",ylab="Catch proportion") {
  txt=txt
  Preds.lower=apply(BootPreds,2,quantile,prob=eps,na.rm=T)
  Preds.upper=apply(BootPreds,2,quantile,prob=1-eps,na.rm=T)
  Pdf=data.frame(len=lenseq,pred=predn,low=Preds.lower,upp=Preds.upper)
  #NB geom_ribbon requires the df to have a y variable, although not used.
  BootGROB=ggplot(data=Pdf,aes(len))+
    #geom_point(data=subset(Tots,subset=c(!is.na(y))),aes(x=lgth,y=y))+
    geom_line(data=Pdf,aes(len,pred))+ylim(0,1)+
    geom_ribbon(data=Pdf,aes(x=len,ymin=low,ymax=upp),alpha=0.2)+
    xlab(xlab)+ylab(ylab)+theme_bw()+
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
  ran.catch=as.data.frame(bind_rows(RanCatchList))
  return(as.data.frame(ran.catch))
}


#' Return a dataframe with permutations of haul treatment
#' @description Randomization of haul treatments (may be more than 2)
#'
#' @param Data Stacked matrix or dataframe of catches in long format. Data are not
#' assumed paired, so could be different length range within each haul.
#' @param haul Name of grouping variable containing the tow/haul/set id number.
#' @param Gear Name of the gear (treatment) indicator variable.
#' @param BlockID Permutation is within blocks
#'
#' @return Dataframe, with randomized gear treatment.
#' @export
#'

RandomizeLong=function(Data,haul="Haul",Gear="Gear",BlockID=NULL) {
  Wk=Data[,c(haul,Gear)]
  if(is.null(BlockID)) Wk[,"BlockID"]="All" else Wk[,"BlockID"]=Data[,BlockID]
  Wk[,haul]=factor(Wk[,haul],levels=unique(Wk[,haul])) #Preserve order
  names(Wk)=c("haul","Gear","BlockID")
  PerSet=Wk %>% group_by(haul,Gear,BlockID) %>% summarize(n=n()) %>% data.frame()
  RanPerSet=PerSet %>% group_by(BlockID) %>%
                 mutate(Gear=sample(Gear)) %>% data.frame()
  RanGear=with(RanPerSet,rep(Gear,n))
  Data[,Gear]=RanGear
  return(Data)
}



#' Return a dataframe with permutations of haul treatment
#' @description Randomization of haul treatments (may be more than 2)
#'
#' @param Data Stacked matrix or dataframe of catches in long format. Data are not
#' assumed paired, so could be different length range within each haul.
#' @param haul Name of grouping variable containing the tow/haul/set id number.
#' @param GearID Name of the gear (treatment) indicator variable.
#'
#' @return Dataframe, with randomized gear treatment.
#' @export
#'

RandomizeLong231105=function(Data,haul="Haul",Gear="Gear") {
  Wk=Data[,c(haul,Gear)]
  Wk[,haul]=factor(Wk[,haul],levels=unique(Wk[,haul])) #Preserve order
  names(Wk)=c("haul","Gear")
  PerSet=Wk %>% group_by(haul,Gear) %>% summarize(n=n()) %>% data.frame()
  RanPerSet=PerSet %>% mutate(Gear=sample(Gear))
  RanGear=with(RanPerSet,rep(Gear,n))
  Data[,Gear]=RanGear
  return(Data)
}



#' Old workhorse
#' @export
#'
Dble.boot230425=function(Data,haul="Haul",BlockID=NULL,
                         Freqs=c("nfine","nwide"),inner.resamp=T,smooth=F) {
  if(smooth) w=c(1,2,1)/4
  Tow=Data[,haul]
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