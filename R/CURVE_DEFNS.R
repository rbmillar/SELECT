## selncurves & propncurves

#' Retention curve definitions.
#' @description The definitions of the retention/selection curves supported by SELECT. The user
#' can edit to add new definitions.
#' @details These curves are used by function propncurves, to calculate the matrix of relative
#' retention probabilities as generated by the R code
#' rmatrix=outer(lens,Meshsize,r,theta).
#' This means that when r is called,
#' both lens and Meshsize will be vectors of length equal to
#' the number of lengthclasses times the number of meshes.
#' @export
#'
selncurves=function(rtype) {
  switch( #Return "logistic" or "richards" for all uses of these curves
    {stype=substr(rtype,4,99)},
   "logistic"={ #Parameters are a and b (i.e., logit link)
      #Geometric similarity is assumed. A control has meshsize of zero.
      function(lens,Meshsize,th) {
      ratio=Meshsize/max(Meshsize) #NB: Baseline is codend
      a=th[1]; b=th[2]
      seln=ifelse(Meshsize==0,1,plogis(a+b*lens/ratio))
      return(seln) } },
   "richards"={ #Parameters are a, b and delta
      #Geometric similarity is assumed. A control has meshsize of zero.
      function(lens,Meshsize,th) {
      ratio=Meshsize/max(Meshsize) #NB: Baseline is codend
      a=th[1]; b=th[2]; delta=exp(th[3])
      seln=ifelse(Meshsize==0,1,plogis(a+b*lens/ratio)^(1/delta))
      return(seln) } },
   "norm.loc"={
       function(lens,Meshsize,th) {
       relsize=Meshsize/min(Meshsize)
       seln=exp(-(lens-th[1]*relsize)^2/(2*th[2]^2))
       return(seln) } },
   "norm.sca"={
       function(lens,Meshsize,th) {
       relsize=Meshsize/min(Meshsize)
       seln=exp(-(lens-th[1]*relsize)^2/(2*th[2]^2*relsize^2))
       return(seln) } },
   "gamma"={
     r=function(lens,Meshsize,th) {
       relsize=Meshsize/min(Meshsize)
       alpha=th[1]; beta=th[2]
       seln=(lens/((alpha-1)*beta*relsize))^(alpha-1)
       seln=seln*exp(alpha-1-lens/(beta*relsize))
       return(seln) } },
   "lognorm"={
       function(lens,Meshsize,th) {
       relsize=Meshsize/min(Meshsize)
       seln=(relsize/lens)*exp(th[1]-th[2]^2/2)
       seln=seln*exp( -(log(lens)-th[1]-log(relsize))^2/(2*th[2]^2) )
       return(seln) } },
   "binorm.sca"={
       function(lens,Meshsize,th) {
       relsize=Meshsize/min(Meshsize)
       seln1=exp(-(lens-th[1]*relsize)^2/(2*th[2]^2*relsize^2))
       seln2=exp(-(lens-th[3]*relsize)^2/(2*th[4]^2*relsize^2))
       p=exp(th[5])/(1+exp(th[5])) #i.e., th[5]=logit(p)
       seln=p*seln1+(1-p)*seln2
       return(seln) } },
   "bilognorm"={
       function(lens,Meshsize,th) {
       relsize=Meshsize/min(Meshsize)
       seln1=(relsize/lens)*exp(th[1]-th[2]^2/2)
       seln1=seln1*exp( -(log(lens)-th[1]-log(relsize))^2/(2*th[2]^2) )
       seln2=(relsize/lens)*exp(th[3]-th[4]^2/2)
       seln2=seln2*exp( -(log(lens)-th[3]-log(relsize))^2/(2*th[4]^2) )
       p=exp(th[5])/(1+exp(th[5])) #i.e., th[5]=logit(p)
       seln=p*seln1+(1-p)*seln2
       return(seln) } },
       stop(paste0('SELECT errror message: ',stype,' not recognised.\n',
        'Possible relative selection types are , "norm.loc", "norm.sca", "gamma",
        "lognorm", "logistic", "richards, "binorm.sca", and "bilognorm". \n',
        'Possible covered-codend and alternative hauls types are "logistic" and "richards" \n'))
   #"binorm.loc"={ #Not yet fully implemented
   #    function(lens,Meshsize,th) {
   #    relsize=Meshsize/min(Meshsize)
   #    seln1=exp(-(lens-th[1]*relsize)^2/(2*th[2]^2))
   #    seln2=exp(-(lens-th[3]*relsize)^2/(2*th[4]^2))
   #    p=exp(th[5])/(1+exp(th[5])) #i.e., th[5]=logit(p)
   #    seln=p*seln1+(1-p)*seln2
   #    return(seln) } },
   #"logistic.L50SR"={ #Parameters are L50 and SR
   #    #Assume geometric similarity when selective cover used
   #    function(lens,Meshsize,th) {
   #    ratio=Meshsize/max(Meshsize) #NB: Baseline is codend
   #    b=2*log(3)/th[2]; a=-b*th[1]
   #    seln=ifelse(Meshsize==0,1,plogis(a+b*lens/ratio))
   #    return(seln) } },
   #"richards.L50SR"={ #Parameters are L50 and SR
   #    #Assume geometric similarity when selective cover used
   #    function(lens,Meshsize,th) {
   #    ratio=Meshsize/max(Meshsize) #NB: Baseline is codend
   #    delta=exp(th[3])
   #    b=(log(0.75^delta/(1-0.75^delta))-log(0.25^delta/(1-0.25^delta)))/th[2]
   #    a=log(0.5^delta/(1-0.5^delta))-b*th[1]
   #    seln=ifelse(Meshsize==0,1,plogis(a+b*lens/ratio)^(1/delta))
   #    return(seln) } },
  )#End of switch
}

#' Relative retention probabilities.
#' @description Returns the function giving the retention probabilities in each of the gears.
#' This takes in to account the selection curve of each gear and the gear configuration
#' (covered-codend, trouser trawl, gillnet).
#' The user can edit to add new definitions.
#' @details These functions are used to calculate the matrix of relative retention probabilities.
#' E.g.,.
#' r=propncurves(rtype)
#' rmatrix=outer(lens,Meshsize,r,theta).
#' The use of the outer function means that when r is called,
#' both lens and Meshsize will be vectors of length equal to
#' the number of lengthclasses times number of meshes.
#' @export
#'
propncurves=function(rtype) {
  switch(
    {dtype=substr(rtype,1,2)},
    "dc"={ selncurves(rtype) },
    "cc"={ #Assume geometric similarity when selective cover used
       function(lens,Meshsize,th) {
        smatrix=matrix(selncurves(rtype)(lens,Meshsize,th),ncol=2)
        seln=c((1-smatrix[,2])*smatrix[,1],smatrix[,2])
        return(seln) } },
    "ec"={ #Assume geometric similarity when selective control used
       function(lens,Meshsize,th) {
        smatrix=matrix(selncurves(rtype)(lens,Meshsize,th),ncol=2)
        psplit=plogis(th[length(th)]) #i.e., th[3]=logit(psplit)
        seln=c((1-psplit)*smatrix[,1],psplit*smatrix[,2])
        return(seln) } },
    "un"={ #Unrestricted parameters of paired hauls
       #The "un" option is for primarily for use with hybrid selectivity
       #and so is not documented or implemented for Estimates()
       function(lens,Meshsize,th) {
        npars=length(th) #Contains pars for both gears, and split
        psplit=plogis(th[npars]) #psplit is last parameter
        nlens=length(lens)
        id1=1:(nlens/2); id2=(nlens/2+1):nlens
        lens1=lens[id1]; lens2=lens[id2]
        Meshsize1=Meshsize[id1]; Meshsize2=Meshsize[id2]
        seln1=selncurves(rtype)(lens1,Meshsize1,th[1:((npars-1)/2)])
        seln2=selncurves(rtype)(lens2,Meshsize2,th[((npars+1)/2):(npars-1)])
        seln=c((1-psplit)*seln1,psplit*seln2)
        return(seln) } },
     stop(paste0('SELECT errror message: ',dtype, ' not recognised.\n',
          'Possible design types are , "dc", "cc", "ec", and "un" \n'))
  )#End of switch
}

