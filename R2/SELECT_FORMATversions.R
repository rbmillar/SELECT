SELECT_FORMAT=function(Df,by=c("haul","lgth"),gear="gear",freq="freq",q.name=NULL,
                       paired=T) {
  UniqueCheck=Df|>group_by(across(all_of(c(by,gear))))|> summarize(m=n(),.groups="keep")
  if(max(UniqueCheck$m)>1)
     stop("\nInput data error: Multiple rows for a unique combination of ",byAll,"\n")
  if(is.null(q.name)) values=freq else values=c(freq,q.name)
  Wk=Df |> select(all_of( c(by,gear,values) ))
  #if(!paired) Wk=Wk |> mutate(uniqRowID=row_number())
  #namePrefix=c("n","q")[1:length(freq)]
  #Wk=Wk |> group_by(across(all_of(by))) |>
   Wk = Wk |>   pivot_wider(names_from=all_of(gear), #names_prefix=namePrefix,
                      values_from=all_of(values), values_fill=0, names_sep="")
  # |> select(-uniqRowID)
  if(!paired) Wk[,gear]=Df[,gear]
  Wk
}
#X=SELECT_FORMAT(Df,by=c("Haul","lgth"),gear="Gear",freq="n",q.name="q"); head(X)

SELECT_FORMAT2=function(Df,by=c("haul","lgth"),gear="gear",freq="freq",q.name=NULL,
                       paired=T) {
  Wk=Df
  if(!paired) Wk=Wk |> mutate(uniqRowID=row_number())
  if(!is.null(q.name)) freq=c(freq,q.name)
  namePrefix=c("n","q")[1:length(freq)]
  Wk=Wk |> group_by(across(all_of(by))) |>
    pivot_wider(names_from=all_of(gear), #names_prefix=namePrefix,
                values_from=all_of(freq), values_fill=0, names_sep="")
  # |> select(-uniqRowID)
  if(!paired) Wk[,gear]=Df[,gear]
  Wk
}

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