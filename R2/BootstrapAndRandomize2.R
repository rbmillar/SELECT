## Bootstrap and randomization functions

#' Bootstrap catch data
#' @description Applies a double bootstrap to data in SELECT format and
#' evaluates the vector valued function `statistic`. The returned value is a 
#' nsim by length(statistic) matrix of bootstrap statistics.
#' @export
bootSELECT2=function(data,var.names,statistic,haul="Haul",nsim=2,
                     block=NULL,gear=NULL,within.resamp=TRUE,...) {
  z=try( statistic(data,...) )
  if(class(z)[1]=="try-error") stop("Error running on actual data")
  #cat("Raw data output:",z,"\n")
  BootMatrix=matrix(NA,nrow=nsim,ncol=length(z))
  Freqs=var.names[-1] #lgth is not needed
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
permSELECT2=function (data,var.names,statistic,haul="Haul", nsim=2,...)
{
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


