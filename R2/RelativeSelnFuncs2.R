#' Fit smooth function to catch share data
#' @description Fit spline to catch share
#'
#' @param Catch Matrix with data in first three columns
#' @param quasi Logical, whether to apply quasibinomial correction
#' @param bs Choice of smoother
#' @param k Dimension of the basis. k=3 in minimum for natural cubic spline.
#' @param m Order of the penalty
#' @param sp Supplied smoothing parameter
#'
#' @return List containing fitted model.
#' @export
fitGAM2=function(var.names,data,scale.names=NULL,useTots=TRUE,
                 bs="cr",k=5,sp=NULL,quasi=FALSE,rm.zeros=TRUE) {
  if(typeof(var.names)!="character")
    stop('SELECT errror: \n Variable names must be character')
  if(!is.null(scale.names) & typeof(scale.names)!="character")
    stop('SELECT errror \n Scaling variable names must be character')
  Data=Raw2Tots(var.names,data,scale.names,useTots)
  vn=var.names
  formla=as.formula( paste0("cbind(",vn[2],",",vn[3],")","~s(",vn[1],",bs=bs,k=k)") )
  if(rm.zeros) Data=Data[Data[,2]+Data[,3]>0,]
  fam=ifelse(quasi,"quasibinomial","binomial")
  GamFit=gam(formla,family=fam,sp=sp,data=Data)
  return(GamFit)
}


Raw2Tots=function(data,var.names,q.names=NULL,scale.names=NULL) {
  if(!is.null(q.names) & !is.null(scale.names)) 
    stop('Raw2Tots error: \n Do not specify both q.names and scale.names ')
  Data=data[,var.names]
  if(!is.null(scale.names)) Data[,-1] = Data[,-1]*data[,scale.names]
  if(!is.null(q.names)) Data[,-1] = Data[,-1]/data[,q.names]
  Data=Data %>% group_by(across(all_of(var.names[1]))) %>%
    summarize(across(all_of(var.names[-1]),sum)) %>% data.frame()
  Data
}



fitGAM2old=function(var.names,data,scale.names=NULL,useTots=T,
                 bs="cr",k=5,sp=NULL,quasi=FALSE,H0=FALSE,rm.zeros=T) {
  if(typeof(var.names)!="character")
    stop('SELECT errror message: \n Variable names must be character')
  if(!is.null(scale.names) & typeof(scale.names)!="character")
    stop('SELECT errror message: \n Scaling variable names must be character')
  Data=data[,var.names]
  if(!is.null(scale.names)) Data[,-1] = Data[,-1]*data[,scale.names]
  if(useTots) Data=Data %>% group_by(across(all_of(var.names[1]))) %>%
    summarize(across(all_of(var.names[-1]),sum)) %>% data.frame()
  formla=as.formula( paste0("y~s(",var.names[1],",bs=bs,k=k)") )
  Data$n=Data[,2]+Data[,3]
  Data$y=Data[,3]/Data$n
  #Data$y[is.na(Data$y)]=0
  if(rm.zeros) Data=Data[Data$n>0,]
  fam=ifelse(quasi,"quasibinomial","binomial")
  if(!H0) GamFit=gam(formla,family=fam,weights=n,sp=sp,data=Data)
  if(H0) {
    GamH0=gam(y~1,family=fam,weights=n,sp=sp,data=Data)
    GamH0.5=gam(y~0,family=fam,weights=n,sp=sp,data=Data)
    GamFit=list(GamH0,GamH0.5)
    }
  return(GamFit)
}

