#' Produce the bootstrap plot with simultaneous confidence intervals
#' @description `BootPlot2` uses ggplot to produce a grob (graphical object)
#' displaying the fitted curve and simultaneous bootstrap confidence intervals.
#' Unlike pointwise intervals, simultaneous intervals ensure that the specified
#' proportion of entire bootstrap curves fall within the bounds.
#'
#' @param BootPreds Matrix with bootstraps by row and fitted values at length in columns, as produced by bootSELECT.
#' @param lenseq Lengths at which fitted values were calculated.
#' @param predn Fitted curve.
#' @param Data If provided the data are added to the plot. The length variable must have the name `lgth` and the catch proportion variable `y`.
#' @param coverage The desired coverage probability for the simultaneous interval. Default is 0.95.
#' @param limits Optional numeric vector of length 2, `c(lower, upper)`, specifying the length range over which simultaneous coverage is computed. Lengths outside this range are excluded when determining whether a bootstrap curve is within bounds. Default is NULL (use all lengths).
#' @param txt Size of text used in the plot axes.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param show.pointwise Logical. If TRUE, also show pointwise intervals as a darker inner band.
#' @return A list with components: `plot` (ggplot GROB), `eps` (the quantile probability used),
#' and `coverage.achieved` (the actual coverage achieved).
#' @details
#' The function finds the value of `eps` such that when bounds are computed as the
#' `eps` and `1-eps` quantiles at each length, approximately `coverage` proportion
#' of bootstrap curves are entirely contained within those bounds.
#'
#' The `limits` argument is useful when extreme lengths have high variability due to
#' sparse data, and you want to prevent these from unduly widening the confidence bands.
#' @export
#'
BootPlot2=function(BootPreds, lenseq, predn, Data=NULL, coverage=0.95, limits=NULL,
                   txt=8, xlab="Length (cm)", ylab="Catch proportion",
                   show.pointwise=FALSE) {

  # Remove rows with NA values
  valid.rows <- complete.cases(BootPreds)
  BootPreds <- BootPreds[valid.rows, , drop=FALSE]
  nsim <- nrow(BootPreds)

  if(nsim < 20) warning("Few valid bootstrap replicates - simultaneous intervals may be unreliable")

  # Determine which columns (lengths) to use for coverage calculation
  if(is.null(limits)) {
    use_cols <- seq_along(lenseq)
  } else {
    use_cols <- which(lenseq >= limits[1] & lenseq <= limits[2])
    if(length(use_cols) == 0) stop("No lengths fall within specified limits")
  }

  # Function to compute coverage for a given eps
  # Coverage = proportion of bootstrap curves entirely within [lower, upper] bounds
  compute_coverage <- function(eps) {
    lower <- apply(BootPreds, 2, quantile, prob=eps, na.rm=TRUE)
    upper <- apply(BootPreds, 2, quantile, prob=1-eps, na.rm=TRUE)
    # Check each bootstrap curve: is it entirely within bounds (for selected columns)?
    within_bounds <- apply(BootPreds, 1, function(curve) {
      all(curve[use_cols] >= lower[use_cols] & curve[use_cols] <= upper[use_cols])
    })
    mean(within_bounds)
  }

  # Binary search to find eps that achieves desired coverage
  # eps=0 gives 100% coverage, eps=(1-coverage)/2 is the pointwise eps
  # Simultaneous eps will always be smaller than pointwise eps
  eps_low <- 0
  eps_high <- (1 - coverage) / 2
  target <- coverage
  tol <- 0.0001

  # Binary search
  max_iter <- 50
  for(i in 1:max_iter) {
    eps_mid <- (eps_low + eps_high) / 2
    cov_mid <- compute_coverage(eps_mid)

    if(abs(cov_mid - target) < tol) {
      break
    }

    # Coverage increases as eps decreases (bounds get wider)
    if(cov_mid < target) {
      # Need wider bounds, decrease eps
      eps_high <- eps_mid
    } else {
      # Bounds too wide, increase eps
      eps_low <- eps_mid
    }
  }

  eps_final <- eps_mid
  coverage_achieved <- cov_mid

  # Compute final bounds
  Preds.lower <- apply(BootPreds, 2, quantile, prob=eps_final, na.rm=TRUE)
  Preds.upper <- apply(BootPreds, 2, quantile, prob=1-eps_final, na.rm=TRUE)

  Pdf <- data.frame(len=lenseq, pred=predn, low=Preds.lower, upp=Preds.upper)

  # Build the plot
  BootGROB <- ggplot(data=Pdf, aes(len)) +
    geom_ribbon(data=Pdf, aes(x=len, ymin=low, ymax=upp), alpha=0.2) +
    geom_line(data=Pdf, aes(len, pred)) + ylim(0,1) +
    xlab(xlab) + ylab(ylab) + theme_bw() +
    theme(axis.text=element_text(size=txt), axis.title=element_text(size=txt)) +
    theme(plot.margin = unit(c(0.75, 0.5, 0.25, 0.5), "cm"))

  # Optionally add pointwise intervals as inner band
  if(show.pointwise) {
    pw_eps <- (1 - coverage) / 2
    Preds.lower.pw <- apply(BootPreds, 2, quantile, prob=pw_eps, na.rm=TRUE)
    Preds.upper.pw <- apply(BootPreds, 2, quantile, prob=1-pw_eps, na.rm=TRUE)
    Pdf$low.pw <- Preds.lower.pw
    Pdf$upp.pw <- Preds.upper.pw
    BootGROB <- BootGROB +
      geom_ribbon(data=Pdf, aes(x=len, ymin=low.pw, ymax=upp.pw), alpha=0.3)
  }

  if(!is.null(Data)) BootGROB <- BootGROB + geom_point(data=Data, aes(x=lgth, y=y))

  # Return list with plot and diagnostics
  list(plot=BootGROB,
       eps=eps_final,
       coverage.achieved=coverage_achieved,
       pointwise.eps=(1-coverage)/2)
}
