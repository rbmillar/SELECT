#' Return a permuted SELECT format dataframe
#' @description Randomization of gear type within each haul.
#' Currently limited to two gears.
#'
#' @param data Matrix or dataframe of catches in SELECT format
#' @param freq.names Character vector giving the names of the two catch frequency variables
#' @param haul Character value giving the haul variable name
#' @param paired Logical. True if the data are paired
#' @param block Character value giving blocking variable. Only used if `paired=FALSE`
#'
#' @return Dataframe, with randomized gear treatment within each haul
#' @export
#'

Randomize=function(data,freq.names=c("n1","n2"),haul="haul",q.names=NULL,
                   paired=TRUE,gear=NULL,block=NULL) {
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
    #Need a numeric haul number since haul order is changed when grouping by block
    Wk$haul=as.numeric(Wk$haul)
    Wk$gear=Wk[,gear]
    if(is.null(block)) Wk$block="All" else Wk$block=data[,block]
    haulgear= Wk %>% group_by(block,haul) %>%
      summarize(haulgear=unique(gear), .groups = "drop_last")
    if(nrow(haulgear)!=nHauls)
      stop("Permutation ERROR: Check data for multiple gear types in a haul")
    #haulgear is grouped by block.
    #Use slice_sample to sample (without replacement) from each block
    #Code works with n=nHauls although it will exceed block size unless block="All"
    permgear = haulgear %>% slice_sample(n=nHauls) %>% pull(haulgear)
    #Now go back to the original haul order
    permuted.gear = permgear[Wk$haul]
    permuted.obs=(Wk$gear!=permuted.gear)
    data[permuted.obs,freq.names]=data[permuted.obs,rev(freq.names)]
    if(!is.null(q.names))
      data[permuted.obs,q.names]=data[permuted.obs,rev(q.names)]
  }
  return(as.data.frame(data))
}
