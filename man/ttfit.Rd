% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrawlFuncs.R
\name{ttfit}
\alias{ttfit}
\title{Fit logistic or Richards selection curve to paired hauls data.}
\usage{
ttfit(
  catch = catchdat,
  type = "logit",
  probs = c(0.25, 0.5, 0.75),
  psplit = NULL,
  x0 = c(-10, 0.4, 0.6),
  delta = 1,
  suff.big = 3,
  nullfit = F,
  plots = T,
  plotfit = T,
  cex = 0.8,
  mkh = 0.07,
  error.bars = FALSE,
  plotlens = NULL,
  details = F,
  xlab = "Length (cm)",
  ylab = "Propn retained",
  verbose = T,
  delta.pen = 0,
  main = c("Proportion of catch in large mesh codend", "Deviance residuals")
)
}
\arguments{
\item{catch}{Matrix with 3 columns. These must be
\enumerate{
\item length (i.e., midpoint of length class);
\item Freq in experimental gear;
\item Freq in nonselective control gear
}}

\item{type}{Character. Use type="rich" to fit Richards curve}

\item{probs}{Numeric to specify retention probabilities for which corresponding
lengths are required}

\item{psplit}{Numeric. If provided the split parameter is fixed at this value}

\item{plotlens}{Numeric. Lengths are which to predict retention probabilities.}

\item{details}{If TRUE, returns additional output}

\item{verbose}{If FALSE then messages are suppressed}
}
\value{
converged Convergence code from nlm minimizer

x Selection curve parameters, a, b, psplit and delta (if Richards curve)

l Log-likelihoods of fitted, null and full models

lens Lengths of 25\%, 50\% and 75\% retention (by default)

sr Selection range

p psplit parameter

CF Overdispersion correction factors, using both deviance and Pearson chisq.
}
\description{
Fit logistic or Richards selection curve to paired hauls data. Pairing would typically be
via twin hauls, alternate hauls or trouser trawls
}
