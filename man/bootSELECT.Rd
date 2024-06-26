% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BOOTSTRAP&PERMUTATION.R
\name{bootSELECT}
\alias{bootSELECT}
\title{Bootstrap catch data}
\usage{
bootSELECT(
  data,
  var.names,
  statistic,
  haul = NULL,
  nsim = 2,
  verbose = T,
  block = NULL,
  gear = NULL,
  within.resamp = TRUE,
  ...
)
}
\description{
Applies a double bootstrap to data in SELECT format and
evaluates the vector valued function \code{statistic}. The returned value is a
nsim by length(statistic) matrix of bootstrap statistics.
}
