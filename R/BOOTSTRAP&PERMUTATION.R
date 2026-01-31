## Bootstrap and randomization functions

#' Bootstrap catch data
#' @description `bootSELECT` applies a hierarchical bootstrap to data
#' in `SELECT` format and evaluates the vector valued function `statistic`.
#' The returned value is a `nsim` by `length(statistic)` matrix of bootstrap statistics.
#' Bootstrapping is first done across `haul`, and then (by default) within `haul`. 
#' 
#' The argument `block` is used to specify bootstrapping across and within blocks. 
#' For example, the blocking variable could be day, or region. See details below.
#'
#' The experiment may be paired of unpaired. In the latter case the `SELECT` format
#' dataframe must contain a variable that indicates the gear type and this variable is
#' used as the `gear=` argument.
#'
#' @param data Stacked matrix or dataframe of catches in SELECT format.
#' @param var.names Character vector of length 3 containing the names of the length variable and catch variables.
#' @param q.names Character vector of length 2 containing the names of the sampling fractions.
#' @param statistic The numeric or vector-valued function to be applied to the bootstrapped data. This function would typically return fit parameters or fitted values.
#' @param haul Name of the grouping variable identifying the haul. This could be a paired-haul in the case of twin or alternate design, in which case both gears must share the same haul identifier.
#' @param nsim Number of bootstrap replicates to be performed.
#' @param paired Logical. This is a required parameter. Set to `TRUE` if the data are paired, `FALSE` otherwise.
#' @param block If specified, name of blocking variable.
#' Bootstrapping is first done over blocks, and then hauls within each block.
#' For example, the blocking variable could be day of deployment.
#' @param gear If specified, name of the gear indicator variable.
#' This is required for use with non-paired data.
#' @param within.resamp Logical. If `TRUE`, then bootstrap resampling is also done at the observatin level within each haul (i.e., double aka, hierarchical bootstrap).
#' @param verbose If set to 0 then messages and the progress bar will be suppressed. If >1 the value of `statistic` for the observed `data` will be printed.
#' @return A matrix of dimension `nsim` by `length(statistic)` containing the bootstrap statistics.
#' @export
bootSELECT=function(data,var.names,q.names=NULL,statistic,haul=NULL,paired=NULL,nsim=2,
                     block=NULL,gear=NULL,within.resamp=TRUE,verbose=1,boot.size=NULL) {
  if(is.null(haul)) stop("The name of the haul variable is required.")
  if(is.null(paired)) stop("The value of paired (TRUE or FALSE) is required")
  if(!paired&is.null(gear))
    stop("\nBOOT ERROR: Analysis of unpaired data requires gear= variable \n")
  if(!paired&!is.null(block))
    cat( crayon::red("WARNING: Bootstrapping over blocks with unpaired data",
        "is inappropriate unless gear efforts are the same in every block\n") )
  statistic.args = list(data=data,var.names=var.names)
  if(!is.null(q.names)) statistic.args$q.names = q.names
  z = try(do.call(statistic, statistic.args))
  if(verbose>1) {
    cat("\n The statistic function applied to the observed data is\n"); print(z) }
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  BootMatrix=matrix(NA,nrow=nsim,ncol=length(z))
  if(verbose) {
    cat(paste("\nStarting a",nsim,"resamples bootstrap...\n"))
    PBar <- txtProgressBar(min = 0, max = nsim, style = 3) }
  for(i in 1:nsim) {
    if(verbose & i%%5==0) setTxtProgressBar(PBar, i)
    bootData=Dble.boot(data,var.names,haul,block,gear,
         paired=paired,within.resamp=within.resamp,boot.size=boot.size)$bootData
    statistic.args$data = bootData
    boot.stat=try(do.call(statistic, statistic.args))
    if(class(boot.stat)[1]=="try-error") {
      cat("\nError running on bootstrap",i,"data\n")
      print(head(bootData)) }
    if(class(boot.stat)[1]!="try-error") BootMatrix[i,]=boot.stat
  }
  if(verbose) close(PBar)
  cat("\nBootstrap successfully completed\n")
  if(any(is.na(BootMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(BootMatrix)
}



#' Return a SELECT format dataframe of hierarchical (double) bootstrapped catch data.
#' @description `Dble.boot` is the double bootstrap function used by `bootSELECT`.
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
#' @return A list object with two components. The first is the bootstrapped dataframe in SELECT format, and the second is a vector containing the bootstrapped `haul` IDs.
#' and
#' @export
#'
Dble.boot=function(data,var.names,haul="Haul",block=NULL,gear=NULL,
                   paired=T,within.resamp=T,smooth=F,boot.size=NULL) {
  Data=as.data.frame(data) #So that code will work with a tibble
  Freqs=var.names[-1] #Catch frequency variables
  if(smooth) w=c(1,2,1)/4
  Sets=Data[,haul]
  uniqSets=unique(Sets)
  #cat("\n",length(uniqSets),"sets to be double bootstrapped")

  #Resample the Sets, subject to gear and block where applicable
  # Bootstrap OVER sets
  if(is.null(block)&is.null(gear)) 
    BootID=if(is.null(boot.size)) SafeSample(uniqSets) else
               SafeSample(uniqSets,size=boot.size)
  # Bootstrap OVER sets WITHIN gear
  if(is.null(block)&!is.null(gear)) {
    GearSetsList=tapply(Sets,Data[,gear],unique)
    GearBootSetsList=if(is.null(boot.size)) lapply(GearSetsList,SafeSample) else
          mapply(SafeSample,GearSetsList,boot.size, SIMPLIFY = FALSE)              
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
  for(j in 1:nRanSets) {
    BootList[[j]] <- Data %>% filter(Sets==BootID[j])
    newHaulID=paste0(j,".",BootList[[j]][,haul])
    BootList[[j]][,haul]=newHaulID }
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
#' @description `bootPlot` uses ggplot to produce a grob (graphical object)
#' displaying the fitted curve and pointwise bootstrap confidence intervals.
#'
#' @param BootPreds Matrix with bootstraps by row and fitted values at length in columns, as produced by bootSELECT.
#' @param lenseqs Lengths at which fitted values were calculated.
#' @param predn Fitted curve
#' @param If provided the data are added to the plot. The length variable must have the name `lgth` and the catch proportion variable `y`.
#' @param eps The quantiles for the lower and upper bootstrap confidence limits are `eps` and 1-`eps`. Default is 95% intervals.
#' @param txt Size of text used in the plot axes.
#' @return ggplot GROB
#' @export
#'
BootPlot=function(BootPreds,lenseq,predn,Data=NULL,coverage=0.95,eps=NULL,txt=8,
                  xlab="Length (cm)",ylab="Catch proportion") {
  if(missing(eps)) eps=(1-coverage)/2
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



#' Produce the bootstrap plot with simultaneous confidence intervals
#' @description `BootPlot2` uses ggplot to produce a grob (graphical object)
#' displaying the fitted curve and simultaneous bootstrap confidence intervals.
#' Unlike pointwise intervals, simultaneous intervals ensure that the specified
#' proportion of entire bootstrap curves fall within the bounds.
#'
#' @param BootPreds Matrix with bootstraps by row and fitted values at length in columns, as produced by bootSELECT.
#' @param lenseq Lengths at which fitted values were calculated.
#' @param predn Fitted curve.
#' @param Data If provided the data are added to the plot. The length variable must have the name `lgth` and the catch proportion variable `y`.
#' @param coverage The desired coverage probability for the simultaneous interval. Default is 0.95.
#' @param limits Optional numeric vector of length 2, `c(lower, upper)`, specifying the length range over which simultaneous coverage is computed. Lengths outside this range are excluded when determining whether a bootstrap curve is within bounds. Default is NULL (use all lengths).
#' @param txt Size of text used in the plot axes.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param show.pointwise Logical. If TRUE, also show pointwise intervals as a darker inner band.
#' @return A list with components: `plot` (ggplot GROB), `eps` (the quantile probability used),
#' and `coverage.achieved` (the actual coverage achieved).
#' @details
#' The function finds the value of `eps` such that when bounds are computed as the
#' `eps` and `1-eps` quantiles at each length, approximately `coverage` proportion
#' of bootstrap curves are entirely contained within those bounds.
#'
#' The `limits` argument is useful when extreme lengths have high variability due to
#' sparse data, and you want to prevent these from unduly widening the confidence bands.
#' @export
#'
BootPlot2=function(BootPreds, lenseq, predn, Data=NULL, coverage=0.95, limits=NULL,
                   txt=8, xlab="Length (cm)", ylab="Catch proportion",
                   show.pointwise=FALSE) {

  # Remove rows with NA values
  valid.rows <- complete.cases(BootPreds)
  BootPreds <- BootPreds[valid.rows, , drop=FALSE]
  nsim <- nrow(BootPreds)

  if(nsim < 99) warning("Few valid bootstrap replicates - intervals may be unreliable")

  # Determine which columns (lengths) to use for coverage calculation
  if(is.null(limits)) {
    use_cols <- seq_along(lenseq)
  } else {
    use_cols <- which(lenseq >= limits[1] & lenseq <= limits[2])
    if(length(use_cols) == 0) stop("No lengths fall within specified limits")
  }

  # Function to compute coverage for a given eps
  # Coverage = proportion of bootstrap curves entirely within [lower, upper] bounds
  compute_coverage <- function(eps) {
    lower <- apply(BootPreds, 2, quantile, prob=eps, na.rm=TRUE)
    upper <- apply(BootPreds, 2, quantile, prob=1-eps, na.rm=TRUE)
    # Check each bootstrap curve: is it entirely within bounds (for selected columns)
    within_bounds <- apply(BootPreds, 1, function(curve) {
      all(curve[use_cols] >= lower[use_cols] & curve[use_cols] <= upper[use_cols])
    })
    mean(within_bounds)
  }
  init.coverage=compute_coverage( (1-coverage)/2 )

  # Binary search to find eps that achieves desired coverage
  # eps=0 gives 100% coverage, eps=(1-coverage)/2 is the pointwise eps
  # It is nimbler to use uniroot (stats) to find eps but it doesn't find exact
  # coverage sometimes, perhaps due to stepwise behaviour of compute_coverage?
  #f=function(eps) compute_coverage(eps)-coverage
  #Target=uniroot(f,lower=0,upper=0.025)
  #eps_final=Target$root
  #coverage_achieved=Target$f.root+coverage
  
  eps_low <- 0
  eps_high <- (1 - coverage) / 2
  target <- coverage
  tol <- 0.0001

  # Binary search
  max_iter <- 50
  for(i in 1:max_iter) {
    eps_mid <- (eps_low + eps_high) / 2
    cov_mid <- compute_coverage(eps_mid)

    if(abs(cov_mid - target) < tol) {
      break
    }

    # Coverage increases as eps decreases (bounds get wider)
    if(cov_mid < target) {
      # Need wider bounds, decrease eps
      eps_high <- eps_mid
    } else {
      # Bounds too wide, increase eps
      eps_low <- eps_mid
    }
  }

  eps_final <- eps_mid
  coverage_achieved <- cov_mid

  # Compute final bounds
  Preds.lower <- apply(BootPreds, 2, quantile, prob=eps_final, na.rm=TRUE)
  Preds.upper <- apply(BootPreds, 2, quantile, prob=1-eps_final, na.rm=TRUE)

  Pdf <- data.frame(len=lenseq, pred=predn, low=Preds.lower, upp=Preds.upper)

  # Build the plot
  BootGROB <- ggplot(data=Pdf, aes(len)) +
    geom_ribbon(data=Pdf[use_cols,], 
                aes(x=len, ymin=low, ymax=upp), alpha=0.2) +
    geom_line(data=Pdf, aes(len, pred)) + ylim(0,1) +
    xlab(xlab) + ylab(ylab) + theme_bw() +
    theme(axis.text=element_text(size=txt), axis.title=element_text(size=txt)) +
    theme(plot.margin = unit(c(0.75, 0.5, 0.25, 0.5), "cm"))

  # Optionally add pointwise intervals as inner band
  if(show.pointwise) {
    pw_eps <- (1 - coverage) / 2
    Preds.lower.pw <- apply(BootPreds, 2, quantile, prob=pw_eps, na.rm=TRUE)
    Preds.upper.pw <- apply(BootPreds, 2, quantile, prob=1-pw_eps, na.rm=TRUE)
    Pdf$low.pw <- Preds.lower.pw
    Pdf$upp.pw <- Preds.upper.pw
    BootGROB <- BootGROB +
      geom_ribbon(data=Pdf, aes(x=len, ymin=low.pw, ymax=upp.pw), alpha=0.3)
  }

  if(!is.null(Data)) BootGROB <- BootGROB + geom_point(data=Data, aes(x=lgth, y=y))

  # Return list with plot and diagnostics
  list(plot=BootGROB,
       initial.pointwise.coverage=coverage,
       initial.simultaneous.coverage=init.coverage,
       adjusted.pointwise.coverage=1-2*eps_final,
       adjusted.simultaneous.coverage=coverage_achieved)
}



#' Sample permutation distribution of a user-defined statistic
#' @description `permSELECT` permutes catch data in SELECT format and returns the value of the user
#' supplied function`statistic`. The returned value is a `nsim` by `length(statistic)` matrix of statistics.
#' 
#' The argument `block` is used to restrict permuting to within blocks. 
#' For example, the blocking variable could be day, or region. See details below.
#'
#' The experiment may be paired of unpaired. In the latter case the `SELECT` format
#' dataframe must contain a variable that indicates the gear type and this variable is
#' used as the `gear=` argument.
#' 
#' @param data Stacked matrix or dataframe of catches in SELECT format.
#' @param var.names Character vector of length 3 containing the names of the length variable and catch variables.
#' #' @param q.names Character vector of length 2 containing the names of the sampling fractions.
#' #' @param statistic The numeric or vector-valued function to be applied to the permuted data. This function would typically be a fit or test statistic.
#' @param haul Name of the grouping variable identifying the haul. This could be a paired-haul in the case of twin or alternate design, in which case both gears must share the same haul identifier.
#' @param nsim Number of permutations  to be performed.
#' @param paired Logical. This is a required parameter. Set to `TRUE` if the data are paired, `FALSE` otherwise.
#' @param block If specified, name of blocking variable. For example, day of deployment. 
#' Permuting is then restricted to being within each block.
#' @param gear If specified, name of the gear indicator variable.
#' This is required for use with non-paired data.
#' @return A matrix of dimension `nsim` by `length(statistic)` containing the bootstrap
#' @export
permSELECT=function (data,var.names,q.names=NULL,statistic,haul="haul",paired=NULL,
                     nsim=2,block=NULL,gear=NULL,verbose=1)
{
  if(is.null(paired)) stop("The value of paired (TRUE or FALSE) is required")
  statistic.args = list(data=data,var.names=var.names)
  if(!is.null(q.names)) statistic.args$q.names = q.names
  z = try(do.call(statistic, statistic.args))
  if(verbose>1) {
    cat("\n The statistic function applied to the observed data is\n"); print(z) }
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  PermMatrix = matrix(NA, nrow = nsim, ncol = length(z))
  if (is.null(haul)) stop("haul is required.")
  if(verbose) {
    cat(paste("\nStarting on", nsim, "permutations...\n"))
    PBar <- txtProgressBar(min = 0, max = nsim, style = 3) }
  for (i in 1:nsim) {
    if (verbose & i%%5 == 0) setTxtProgressBar(PBar, i)
    permData = Randomize(data,var.names,q.names,haul=haul,
                         paired=paired,gear=gear,block=block)
    statistic.args$data = permData
    perm.stat = try(do.call(statistic, statistic.args))
    if (class(perm.stat)[1] == "try-error") {
      cat("\nError running on permutation", i, "data\n")
      print(head(permData))
    }
    if (class(perm.stat)[1] != "try-error")
      PermMatrix[i, ] = perm.stat
  }
  if(verbose) close(PBar)
  cat("\nPermutations successfully completed\n")
  if (any(is.na(PermMatrix)))
    cat("CAUTION: Some fits did not converge - please check for NAs in output.\n")
  invisible(PermMatrix)
}



#' Returns a permuted SELECT format dataframe
#' @description `Randomize` randomly permutes gear type. If the experiment is paired haul then permutation is solely within each paired haul. If unpaired then permutation is across all hauls, possibly within blocks. Currently limited to two gears.
#'
#' @param data Matrix or dataframe of catches in SELECT format
#' @param freq.names Character vector giving the names of the two catch frequency variables
#' @param haul Character value giving the haul variable name. This must be unique for all hauls.
#' @param paired Logical. True if the data are paired
#' @param block Character value giving blocking variable. Only used if `paired=FALSE`
#'
#' @return Dataframe, with randomized gear treatment within each haul
#' @export
#'
Randomize=function(data,var.names,q.names=NULL,haul="haul",
                   paired=TRUE,gear=NULL,block=NULL) {
  if(paired==T&(!is.null(block)|!is.null(gear)))  cat("\n NOTE: gear and block
        variables will be ignored since permutation is within gear pairs\n")
  if(paired==F&is.null(gear)) Stop("\n ERROR: gear name is required
        since data are unpaired and each row must be identified by gear type \n")
  Wk=data.frame(data)
  Wk$haul=as.factor(Wk[,haul])
  nHauls=length(unique(Wk$haul))
  freq.names=var.names[-1]
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
    #Identify the permuted rows in the data frame, p.obs
    p.obs=(Wk$haul %in% permuted.hauls)
    #Create new variables Permuted (T or F) and permGear (gear name)
    data$Permuted=p.obs
    Gears=unique(Wk$gear)
    permGear=paste0("permuted.",gear)
    data[,permGear]=data[,gear]
    data[p.obs,permGear]=ifelse(data[p.obs,permGear]==Gears[1],Gears[2],Gears[1])
    #Permute the frequencies and sampling fractions
    data[p.obs,freq.names]=data[p.obs,rev(freq.names)]
    if(!is.null(q.names))
      data[p.obs,q.names]=data[p.obs,rev(q.names)]
  }
  return(as.data.frame(data))
}



#' Compute statistics suitable for permutation testing
#' @description `SplineStatistics` calculates a variety of test statistics for use
#' with permutation testing. The appropriate test statistic(s) to use will depend on
#' the questions of interest. It does NOT assume that the SplineSELECT fit
#' returns a log-likelihood and so still work when a quasi=T fit is used.
#' @details
#' This function returns several statistics, including:
#' \itemize{
#'   \item `DevExpl` - the proportion of deviance explained by the model
#'   \item `EqualDevExpl` - the proportion of deviance explained by the model,
#'   but with null model fitted with a constant catch share of `Equal`
#'   \item `null` - the log-likelihood of the null model
#'   \item `Equal` - the log-likelihood of the model with a constant catch share of `Equal`
#'   \item `full` - the log-likelihood of the full model
#'   \item `model` - the log-likelihood of the fitted model
#'   \item `LRT` - the likelihood ratio test statistic comparing the fitted spline model
#'    to the null model. Used to test for evidence of a length effect.
#'   \item `EqualLRT` - the likelihood ratio test statistic comparing the fitted model
#'   to the model with a constant catch share of `Equal`
#'   \item `RatioPropnMLS` - the ratio (across the two gears) of the proportion of catch
#'   that is of MLS or bigger.
#'   \item `PropnGear2` - the average catch share in the second gear
#' }
#'
#' @export
SplineStatistics=function(SplineFit,Equal=0.5,MLS=NULL) {
  if(!"gam" %in% class(SplineFit)) stop("Fitted model must be a gam")
  D=summary(SplineFit)$dev.expl
  yhat=fitted(SplineFit)
  nBA=SplineFit$model[[1]] #Numbers in gear B (column 1) and A (column 2)
  ybar=sum(nBA[,1])/sum(nBA); #verage catch share
  n=nBA[,1]+nBA[,2]
  y=nBA[,1]
  RatioPropnMLS=NA
  if(!is.null(MLS)){
    CL=SplineFit$model[[2]]
    TotsBA=apply(nBA,2,sum)
    MLSpropns=apply(nBA[CL>=MLS,],2,sum)/TotsBA
    RatioPropnMLS=MLSpropns[1]/MLSpropns[2]
    names(RatioPropnMLS)=NULL }
  EqualLLhood=sum( dBinom(y,n,Equal,log=T) )
  NullLLhood=sum( dBinom(y,n,ybar,log=T) )
  SplineLLhood=sum( dBinom(y,n,yhat,log=T) )
  FullLLhood=sum( dBinom(y,n,ifelse(n>0,y/n,0),log=T) )
  LRT=2*(SplineLLhood-NullLLhood)
  EqualLRT=2*(SplineLLhood-EqualLLhood)
  ModelDev=2*(FullLLhood-SplineLLhood)
  NullDev=2*(FullLLhood-NullLLhood)
  EqualDev=2*(FullLLhood-EqualLLhood)
  D=1-ModelDev/NullDev #Or, D=summary(SplineFit)$dev.expl
  EqualD=1-ModelDev/EqualDev
  Stats=c(DevExpl=D,EqualDevExpl=EqualD,null=NullLLhood, Equal=EqualLLhood,
          full=FullLLhood,model=SplineLLhood,LRT=LRT,EqualLRT=EqualLRT,
          RatioPropnMLS=RatioPropnMLS,PropnGear2=ybar)
}




#' Calculate the permutation p-value
#'
#' @description Convenience function to evaluate the permutation p-value by comparing the observed value of the statistic, `ObsStat`, to the permuted values in `PermOut`.
#' @export
#'
permPval=function(ObsStat,PermOut,signif="greater",includeObs=T,eps=1e-10) {
  nsim=length(PermOut)
  permDiff=PermOut-ObsStat+ifelse(signif=="greater",eps,-eps)
  permSum=ifelse(signif=="greater",sum(permDiff>=0),sum(permDiff<=0))
  if(!includeObs) permPval=permSum/nsim else
    permPval=(permSum+1)/(nsim+1)
  return(permPval)
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
SafeSample=function(x,size=length(x),replace=T) {
  if(length(x)>1) x=sample(x,size,replace=replace)
  return(x)
}
