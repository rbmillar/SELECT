
#######################
## Deprecated functions
#######################


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


#' Return a permuted SELECT format dataframe
#' @description Randomization of gear type within each haul
#'
#' @param catch Stacked matrix or dataframe of catches in SELECT format
#' @param varnames Length-3 character vector containing name of haulID variable, and 2 catch freq vars
#'
#' @return Dataframe, with randomized gear treatment within each haul
#' @export
#'

Randomize240613=function(catch,varnames=c("Haul","nwide","nfine")) {
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
