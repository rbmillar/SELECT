% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HELPER_FUNCTIONS.R
\name{ModelCheck}
\alias{ModelCheck}
\title{Provide a diagnostic summary of the SELECT model fit}
\usage{
ModelCheck(
  fit,
  minE = 0,
  xlab = "Length (cm)",
  ylab = "Propn in exptl gear",
  print.out = T,
  plots = c(T, T),
  plotlens = NULL,
  pex = 1,
  ...
)
}
\description{
Provides log-likelihoods, deviances,
over-dispersion correction factors (if minE>0), and a plot of deviance residuals.
}
