% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xbadzupa.R
\name{xbadzupa}
\alias{xbadzupa}
\title{Perform XBADZUPA analysis}
\usage{
xbadzupa(
  dat,
  Nbootstrap = 400,
  ci = 0.9,
  showDetectedPeaks = TRUE,
  ignoreWarn = TRUE
)
}
\arguments{
\item{dat}{The data for peak detection.}

\item{Nbootstrap}{The number of bootstrap samples to generate.}

\item{ci}{The confidence level for peak detection.}

\item{showDetectedPeaks}{Whether to display detected peaks.}

\item{ignoreWarn}{Whether to ignore warning messages.}
}
\description{
This function performs XBADZUPA peak detection analysis on the provided data.
}
\examples{
set.seed(1234); dat <- c(rnorm(200, 4), rnorm(50, 8, 2))
xbadzupa(dat, 100)

}
