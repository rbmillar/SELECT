
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

SELECT_FORMAT=function(Df,by=c("TowID","lgth"),gear="gear",freq="freq") {
  wk=split(Df,Df[,gear])
  ngear=length(wk)
  nby=length(by)
  freq.names=paste0("n",names(wk))
  Stacked.df=full_join(wk[[1]][,c(by,freq)],wk[[2]][,c(by,freq)],by=by)
  if(ngear>2) {
    for(k in 3:ngear)
      Stacked.df=full_join(Stacked.df,wk[[k]][,c(by,freq)],by=by) }
  names(Stacked.df)[nby+(1:ngear)]=freq.names
  #Reshuffle columns to length is in column 1
  Stacked.df=Stacked.df[,c(nby+(0:ngear),1:(nby-1))]
  Stacked.df[is.na(Stacked.df)] <- 0
  Stacked.df
}
