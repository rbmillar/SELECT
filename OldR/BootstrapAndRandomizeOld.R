#Double bootstrap of multihaul data in SELECT format

#Crude weighted average (0.25,0.5,0.25) with immediate neighbours
WgtAvg=function(y,w=c(0.25,0.5,0.25)) {
  n=length(y)
  y.left=c(y[1],y)
  y.right=c(y,y[n])
  wgt.y=w[2]*y+w[1]*y.left[1:n]+w[3]*y.right[2:(n+1)]
  wgt.y
}
#' Double bootstrap function used by SELECT package
#'
#' @param Data Stacked matrix or dataframe of catches in SELECT format,
#' with lengthclass in first column. Remaining columns are raw catch frequencies
#' @param Set.id Grouping variable containing the tow/haul/set id number.
#' Vector of same length as number of rows in Data.
#' @param smooth Smooth at within-haul phase to avoid losing degrees of freedom (from increasing number of zero freqs)
#'
#' @return Dataframe of double-bootstrapped freqs
#' @export
#'

Dble.boot=function(Data,SetID="Haul",Freqs=c("nfine","nwide"),smooth=T) {
  Tow=Data[,SetID]
  uniqTows=unique(Tow)
  nTows=length(uniqTows)
  #cat("\n",nTows,"hauls to be double bootstrapped")
  BootCatchList=as.list(1:nTows)
  BootID=sample(uniqTows,nTows,replace=T)
  for(j in 1:nTows) { BootCatchList[[j]] <- Data %>% filter(Tow==BootID[j]) }
  bootData=as.data.frame(rbindlist(BootCatchList))
  if(smooth) w=c(1,2,1)/4 else w=c(0,1,0)
  #Parametric bootstrap within tows
  nObs=nrow(bootData)
  for(j in 1:length(Freqs))
    bootData[,Freqs[j]]=rpois(nObs,WgtAvg(bootData[,Freqs[j]],w=w))
  return(list(bootData=bootData,BootID=BootID))
}


Dble.boot2=function(Data,SetID,smooth=T) {
  uniqTows=unique(SetID)
  nTows=length(uniqTows)
  #cat("\n",nTows,"hauls to be double bootstrapped")
  BootDataList=as.list(1:nTows)
  BootID=sample(uniqTows,nTows,replace=T)
  for(j in 1:nTows) { BootDataList[[j]] <- Data %>% filter(SetID==BootID[j]) }
  bootData=as.data.frame(rbindlist(BootDataList))
  if(smooth) w=c(1,2,1)/4 else w=c(0,1,0)
  #Parametric bootstrap within tows
  nObs=nrow(bootData)
  #Assumes only two tows for now, with freqs in columns 2&3
  for(j in 2:3) {bootData[,j]=rpois(nObs,WgtAvg(bootData[,j],w=w))}
  return(list(bootData=bootData,BootID=BootID))
}

#' Double bootstrap function
#'
#' @param catch Stacked matrix or dataframe of catches in SELECT format
#' @param varnames Length-3 vector containing name of haulID variable, and 2 catch freq vars
#' @param smooth Smooth at within-haul phase to avoid losing degrees of freedom (from increasing number of zero freqs)
#'
#' @return Dataframe of double-bootstrapped freqs
#' @export
#'

Dble.bootstrap=function(catch,varnames=c("Haul","nwide","nfine"),smooth=T) {
  Tow=catch[,varnames[1]]
  uniqTows=unique(Tow)
  nTows=length(uniqTows)
  #cat("\n",nTows,"hauls to be double bootstrapped")
  BootCatchList=as.list(1:nTows)
  BootID=sample(uniqTows,nTows,replace=T)
  for(j in 1:nTows) {
    BootCatchList[[j]] <- catch %>% filter(Tow==BootID[j])
    BootCatchList[[j]]$Tow <- j
    BootCatchList[[j]]$Orig.Tow <- BootID[j]
  }
  boot.catch=as.data.frame(rbindlist(BootCatchList))
  if(smooth) w=c(1,2,1)/4 else w=c(0,1,0)
  #Parametric bootstrap within tows
  nObs=nrow(boot.catch)
  boot.catch[,varnames[2]]=rpois(nObs,WgtAvg(boot.catch[,varnames[2]],w=w))
  boot.catch[,varnames[3]]=rpois(nObs,WgtAvg(boot.catch[,varnames[3]],w=w))
  return(as.data.frame(boot.catch))
}

#' Randomization of gear type within each haul
#'
#' @param Stacked matrix or dataframe of catches in SELECT format
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