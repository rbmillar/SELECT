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
  paired = NULL,
  nsim = 2,
  block = NULL,
  gear = NULL,
  within.resamp = TRUE,
  verbose = 1,
  ...
)
}
\arguments{
\item{data}{Stacked matrix or dataframe of catches in SELECT format.}

\item{var.names}{Character vector of length 3 containing the names of the length variable and catch variables.}

\item{statistic}{The numeric or vector-valued function to be applied to the bootstrapped data. This function would typically return fit parameters or fitted values.}

\item{haul}{Name of the grouping variable identifying the haul. This could be a paired-haul in the case of twin or alternate design, in which case both gears must share the same haul identifier.}

\item{paired}{Logical. This is a required parameter. Set to \code{TRUE} if the data are paired, \code{FALSE} otherwise.}

\item{nsim}{Number of bootstrap replicates to be performed.}

\item{block}{If specified, name of blocking variable.
Bootstrapping is first done over blocks, and then hauls within each block.
For example, day of deployment.}

\item{gear}{If specified, name of the gear indicator variable.
This is required for use with non-paired data.}

\item{within.resamp}{Logical. If \code{TRUE}, then bootstrap resampling is also done at the observatin level within each haul (i.e., double aka, hierarchical bootstrap).}

\item{verbose}{If set to 0 then messages and the progress bar will be suppressed. If >1 the value of \code{statistic} for the observed \code{data} will be printed.}

\item{...}{Other parameters to be passed to the \code{statistic} function. E.g.,
q.names (sampling fractions) if \code{statistic} has been defined to take \code{q.names} as an argument.}
}
\value{
A matrix of dimension \code{nsim} by \code{length(statistic)} containing the bootstrap statistics.
}
\description{
Applies a hierarchical (double) bootstrap to data in SELECT format and
evaluates the vector valued function \code{statistic}. The returned value is a
\code{nsim} by \code{length(statistic)} matrix of bootstrap statistics.
}
