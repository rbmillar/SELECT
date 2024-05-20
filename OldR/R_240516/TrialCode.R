#To demonstrate use of model formula rather

Poly2=function(fo,total,data,Quasi=F,wgt="AICc",All=TRUE) {
  if(wgt!="AICc") wgt="AIC" #Only two possibilities for now
  if(All) Sub=expression(eval(TRUE))
  if(!All) Sub=expression( dc(I(x^0),x,I(x^2),I(x^3),I(x^4)) )
  yvar=deparse(fo[[2]])
  xvar=deparse(fo[[3]])
  totvar=as.character(substitute(total))
  cat(yvar,xvar,totvar,"\n")
  cat(typeof(yvar),typeof(xvar),typeof(totvar))
  Catch=data[,c(xvar,yvar,totvar)]
  Catch
}

xxx=Poly2(y~lgth,n,TotER.df)
