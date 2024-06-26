% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BOOTSTRAP&PERMUTATION.R
\name{BootPlot}
\alias{BootPlot}
\title{Produce the bootstrap plot}
\usage{
BootPlot(
  BootPreds,
  lenseq,
  predn,
  Data = NULL,
  eps = 0.025,
  txt = 8,
  xlab = "Length (cm)",
  ylab = "Catch proportion"
)
}
\arguments{
\item{BootPreds}{Matrix with bootstrap by row and fitted values at length in columns.}

\item{predn}{Fitted curve}

\item{lenseqs}{Lengths at which fitted values are calculated.}
}
\value{
ggplot GROB
}
\description{
ggplot showing fitted curve and bootstrap bounds
}
