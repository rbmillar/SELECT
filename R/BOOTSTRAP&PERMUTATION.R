## Bootstrap and randomization functions

#' Bootstrap catch data
#' @description Applies a double bootstrap to data in SELECT format and
#' evaluates the vector valued function `statistic`. The returned value is a
#' nsim by length(statistic) matrix of bootstrap statistics.
#'
#' @param data Stacked matrix or dataframe of catches in SELECT format,
#' @param haul Name of grouping variable containing the tow/haul/set id number.
#' Vector of same length as number of rows in Data.
#' @param block If specified, name of blocking variable.
#' Bootstrapping is first done over blocks, and then sets within each block.
#' @param gear If specified, name of the gear indicator variable.
#' This restricts bootstrapping of sets to within each gear and is
#' for use with non-paired data.
#' @param paired Logical True if the data are paired
#' @param verbose If set to 0 will suppress messages and progress bar.
#' If >1 will also print the value of `statistic` for the observed `data`.
#' @param ... Other parameters to be passed to the `statistic` function. E.g.,
#' q.names (sampling fractions) if applicable.
#' @export
bootSELECT=function(data,var.names,statistic,haul=NULL,nsim=2,paired=NULL,
                     block=NULL,gear=NULL,within.resamp=TRUE,verbose=1,...) {
  if(is.null(haul)) stop("The name of the haul variable is required.")
  if(is.null(paired)) stop("The value of paired (TRUE or FALSE) is required")
  if(!paired&is.null(gear))
    stop("\nBOOT ERROR: Analysis of unpaired data requires gear= variable \n")
  if(!paired&!is.null(block))
    cat( crayon::red("WARNING: Bootstrapping over blocks with unpaired data",
        "is inappropriate unless gear efforts are the same in every block\n") )
  #data=obj$data; var.names=obj$var.names; z=try( statistic(data,...) )
  z=try( statistic(data,var.names,...) )
  if(verbose>1) {
    cat("\n The statistic function applied to the observed data is\n"); print(z) }
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  BootMatrix=matrix(NA,nrow=nsim,ncol=length(z))
  Freqs=var.names[-1] #lgth is not needed
  if(verbose) {
    cat(paste("\nStarting a",nsim,"resamples bootstrap...\n"))
    PBar <- txtProgressBar(min = 0, max = nsim, style = 3) }
  for(i in 1:nsim) {
    if(verbose & i%%5==0) setTxtProgressBar(PBar, i)
    bootData=Dble.boot(data,haul,block,gear,Freqs,
                       paired=paired,within.resamp)$bootData
    boot.stat=try( statistic(bootData,var.names,...) )
    if(class(boot.stat)[1]=="try-error") {
      cat("\nError running on bootstrap",i,"data\n")
      print(head(bootData)) }
    if(class(boot.stat)[1]!="try-error") BootMatrix[i,]=boot.stat
  }
  if(verbose) {
    close(PBar)
    cat("\nBootstrap successfully completed\n") }
  if(any(is.na(BootMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(BootMatrix)
}



#' Permute catch data
#' @description Permutes data in SELECT format and returns the value of the user
#' supplied function`statistic`in a nsim by length(statistic) matrix
#' @export
permSELECT=function (data,var.names,statistic,haul="haul",paired=NULL, nsim=2,
                     gear=NULL,block=NULL,...)
{
  #data=obj$data; var.names=obj$var.names; z=try( statistic(data,...) )
  if(is.null(paired)) stop("The value of paired (TRUE or FALSE) is required")
  z = try(statistic(data,var.names,...))
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
    permData = Randomize(data,freq.names=Freqs,haul=haul,
                         paired=paired,gear=gear,block=block,...)
    perm.stat = try(statistic(permData,var.names,...))
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
#' @description Double bootstrap function used by bootSELECT function
#'
#' @param data Stacked matrix or dataframe of catches in SELECT format,
#' with lengthclass in first column. Remaining columns are raw catch frequencies
#' @param haul Name of grouping variable containing the tow/haul/set id number.
#' Vector of same length as number of rows in Data.
#' @param block If specified, name of blocking variable.
#' Bootstrapping is first done over blocks, and then sets within each block.
#' @param gear If specified, name of the gear indicator variable.
#' This restricts bootstrapping of sets to within each gear and is intended to be used
#' for non-paired data.
#' @param paired Logical. True if the data are paired
#' @param within.resamp If F, no resampling at observation level.
#' @param smooth Smooth at within-haul phase to avoid losing degrees of freedom (from increasing number of zero freqs)
#'
#' @return Dataframe of double-bootstrapped freqs
#' @export
#'

Dble.boot=function(data,haul="Haul",block=NULL,gear=NULL,
                   Freqs=c("nfine","nwide"),paired=T,within.resamp=T,smooth=F) {
  Data=as.data.frame(data) #So that code will work with a tibble
  if(smooth) w=c(1,2,1)/4
  Sets=Data[,haul]
  uniqSets=unique(Sets)
  nSets=length(uniqSets)
  #cat("\n",nSets,"sets to be double bootstrapped")
  RanList=list(uniqSets)

  #Resample the Sets, subject to gear and block where applicable
  # Bootstrap OVER sets
  if(is.null(block)&is.null(gear)) BootID=sample(uniqSets,nSets,replace=T)
  # Bootstrap OVER sets WITHIN gear
  if(is.null(block)&!is.null(gear)) {
    GearSetsList=tapply(Sets,Data[,gear],unique)
    GearBootSetsList=lapply(GearSetsList,SafeSample,replace=T)
    BootID=unlist(GearBootSetsList)  }
  #Bootstrap OVER blocks, then OVER sets or sets WITHIN gear
  if(!is.null(block)) {
    BlockSetsList=tapply(Sets,Data[,block],unique)
    nBlocks=length(BlockSetsList)
    BootBlockSetsList=BlockSetsList[sample(1:nBlocks,replace=T)]
    # Sample the Sets within block
    if(is.null(gear)) {
      BootBlockBootSetList=lapply(BootBlockSetsList,SafeSample,replace=T)
      BootID=unlist(BootBlockBootSetList) }
    if(!is.null(gear)) {
      GearSetsList=tapply(Sets,Data[,gear],unique)
      nGears=length(GearSetsList)
      WkList=as.list(1:nGears)
      for(j in 1:nGears) {
        WkList[[j]]=lapply(BootBlockSetsList,function(x) x[x %in% GearSetsList[[j]] ])
        WkList[[j]]=lapply(WkList[[j]],SafeSample,replace=T) }
      BootID=unlist(WkList) }
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
#' @description Randomization of gear type within each haul.
#' Currently limited to two gears.
#'
#' @param data Matrix or dataframe of catches in SELECT format
#' @param freq.names Character vector giving the names of the two catch frequency variables
#' @param haul Character value giving the haul variable name. This must be UNIQUE for all hauls.
#' @param paired Logical. True if the data are paired
#' @param block Character value giving blocking variable. Only used if `paired=FALSE`
#'
#' @return Dataframe, with randomized gear treatment within each haul
#' @export
#'

Randomize=function(data,freq.names=c("n1","n2"),haul="haul",
                   paired=TRUE,gear=NULL,block=NULL,q.names=NULL) {
  if(paired==T&(!is.null(block)|!is.null(gear)))  cat("\n NOTE: gear and block
        variables will be ignored since permutation is within gear pairs\n")
  if(paired==F&is.null(gear)) Stop("\n ERROR: gear name is required
        since data are unpaired and each row must be identified by gear type \n")
  Wk=data.frame(data)
  Wk$haul=as.factor(Wk[,haul])
  nHauls=length(unique(Wk$haul))
  if(paired) {
    uniqHauls=unique(Wk[,haul])
    permList=as.list(1:nHauls)
    for(j in 1:nHauls) {
      permList[[j]] = Wk %>% filter(haul==uniqHauls[j])
      jPerm=sample(1:2)
      permList[[j]][,freq.names] <- permList[[j]][,freq.names[jPerm]]
      if(!is.null(q.names))
        permList[[j]][,q.names] <- permList[[j]][,q.names[jPerm]]
    }
  data=bind_rows(permList)
  }
  if(!paired) {
    Wk$gear=Wk[,gear]
    if(is.null(block)) Wk$block="All" else Wk$block=data[,block]
    haulgrp= Wk %>% group_by(block,haul) %>%
      summarize(grp=unique(gear), .groups = "drop_last")
    if(nrow(haulgrp)!=nHauls)
      stop("Permutation ERROR: Check data for multiple gear types in a haul")
    #haulgear tibble is still grouped by block.
    #Use slice_sample to sample (without replacement) from each block
    #Code works with n=nHauls even though it exceeds block size unless block="All"
    permgrp = haulgrp %>% slice_sample(n=nHauls)
    #Identify which hauls are permuted
    permuted.hauls=haulgrp$haul[haulgrp$grp!=permgrp$grp]
    #Identify the permuted rows in the data frame, and permute
    permuted.obs=(Wk$haul %in% permuted.hauls)
    data$Permuted=permuted.obs
    data[permuted.obs,freq.names]=data[permuted.obs,rev(freq.names)]
    if(!is.null(q.names))
      data[permuted.obs,q.names]=data[permuted.obs,rev(q.names)]
  }
  return(as.data.frame(data))
}

#' @export
#'
permPval=function(ObsStat,PermOut,signif="greater",includeObs=FALSE) {
  nsim=length(PermOut)
  permDiff=PermOut-ObsStat
  permSum=ifelse(signif=="greater",sum(permDiff>=0),sum(permDiff<=0))
  if(!includeObs) permPval=permSum/nsim else
    permPval=(permSum+1)/(nsim+1)
  return(permPval)
}

