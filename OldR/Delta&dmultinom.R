#Modification of deltamethod (msm) that also returns function value
delta.method=function (g, mean, cov, ses = TRUE)
{
  cov <- as.matrix(cov)
  n <- length(mean)
  if (!is.list(g))
    g <- list(g)
  if ((dim(cov)[1] != n) || (dim(cov)[2] != n))
    stop(paste("Covariances should be a ", n, " by ",
               n, " matrix"))
  syms <- paste("x", 1:n, sep = "")
  for (i in 1:n) assign(syms[i], mean[i])
  gval = t(sapply(g, function(form) { as.numeric(eval(deriv(form, syms))) } ))
  gdashmu <- t(sapply(g, function(form) {
    as.numeric(attr(eval(deriv(form, syms)), "gradient")) }))
  new.covar <- gdashmu %*% cov %*% t(gdashmu)
  if (ses) { list(value=as.vector(gval),se=sqrt(diag(new.covar))) }
  else list(value=as.vector(gval),covar=new.covar)
}


#' @export
dmultinomial=function (x, size = NULL, prob, log = FALSE)
{
  if (length(x) == 0)
    return(numeric(0))
  if (is.vector(prob))
    prob <- t(prob)
  if (is.vector(x))
    x <- t(x)
  K <- ncol(prob)
  maxn <- max(nx <- nrow(x), np <- nrow(prob))
  if (ncol(x) != K)
    stop("x[] and prob[] must be equal length vectors or equal col matrix.")
  if (nx != maxn)
    x <- matrix(t(x), ncol = K, nrow = maxn, byrow = TRUE)
  if (np != maxn)
    prob <- matrix(t(prob), ncol = K, nrow = maxn, byrow = TRUE)
  N <- apply(x, 1, sum)
  if (is.null(size))
    size <- N
  size <- rep(size, length = maxn)
  if (any(size != N))
    stop("at least one size != sum(x), i.e. one is wrong")
  if (any(prob < 0) || any((s <- apply(prob, 1, sum)) == 0))
    stop("probabilities cannot be negative nor all 0.")
  prob <- prob/s
  if (any(as.integer(x + 0.5) < 0))
    stop("'x' must be non-negative")
  r <- sapply(1:maxn, function(y) {
    xx <- as.integer(x[y, ] + 0.5)
    pp <- prob[y, ]
    i0 <- pp == 0
    if (any(i0)) {
      if (any(xx[i0] != 0))
        return(0)
      xx <- xx[!i0]
      pp <- pp[!i0]
    }
    return(lgamma(size[y] + 1) + sum(xx * log(pp) - lgamma(xx + 1)))
  })
  if (log)
    return(r)
  else return(exp(r))
}


#' @export
byrow_dmultinom=function(x, prob, log = FALSE) {
  x=as.matrix(x); prob=as.matrix(prob)
  if(nrow(x)!=nrow(prob)|ncol(x)!=ncol(prob))
    stop("x[] and prob[] must have same dimension.")
  if(ncol(x)==1) stop("x[] must have at least 2 columns")
    k=ncol(x)
    rowN=apply(x,1,sum)
    r=sapply(1:nrow(x), function(i) {
      xx <- x[i,]
      pp <- prob[i, ]
      i0 <- pp == 0
      if (any(i0)) {
        if (any(xx[i0] != 0))
          return(0)
        xx <- xx[!i0]
        pp <- pp[!i0]
      }
      return(lgamma(rowN[i] + 1) + sum(xx*log(pp) - lgamma(xx+1)))
    })
    if (log)
      return(r)
    else return(exp(r))
}
#fullfit.l=sum(byrow_dmultinom(Counts,prob=CountPropns,log=TRUE))


