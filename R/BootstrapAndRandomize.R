## Bootstrap and randomization functions

#' Implement bootstrap
#' @description Apply the double bootstrap and return a statistic from each
#' bootstrap sample in a matrix.
#' @export
boot.SELECT=function(Data,statistic,N,SetID="Haul",Freqs=c("nfine","nwide"),...) {
  z=try( statistic(Data,...) )
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  #cat("Raw data output:",z,"\n")
  BootMatrix=matrix(NA,nrow=N,ncol=length(z))
  cat(paste("Bootstrap data: Starting a",N,"resamples bootstrap...\n"))
  if(is.null(SetID)) stop("SetID (set id variable) is required.")
  PBar <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {
    if(i%%5==0) setTxtProgressBar(PBar, i)
    bootData=Dble.boot(Data,SetID,Freqs,smooth=F)$bootData
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



#' Permutation test for no difference in selectivity of two gear
#' @description Randomization of gear type within each haul
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
