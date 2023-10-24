#' Get the best algorithm for density estimation using cross-validation.
#'
#' This function uses cross-validation to determine the best density estimation
#' algorithm among the provided options.
#'
#' @param d The data for density estimation.
#' @param algorithms A list of density estimation functions to consider.
#' @return function, the best density estimation algorithm.
#' @examples
#' set.seed(1234); dat <- c(rnorm(200, 4), rnorm(50, 8, 2))
#' bestAlgorithm <- getBestAlgorithmByCv(dat)
#' plot(bestAlgorithm(dat))
#' 
#' @export 
getBestAlgorithmByCv <- function(d, algorithms = list(
      densityADEBA = function(x) density.adeba(x, out=0.00, n=512),
      densityBOTEV = function(x) IsoplotR::kde(x, plot=FALSE)
  )) {

    cv <- Cv$new(algorithms, d)
    cv$usingLogStar = TRUE
    res <- cv$exeCVinAllAlgorithm()

    algorithm_score <- res$algorithm %>% 
    unique %>% 
    sapply(function (x) {
      mean(res[res$algorithm==x,]$score)
    })
    return(algorithms[[names(which.min(algorithm_score))]])
}

#' Get bootstrap results for a density estimation algorithm.
#'
#' This function performs bootstrapping to obtain confidence intervals and other
#' information for a density estimation algorithm.
#'
#' @param dat The data for density estimation.
#' @param densFn The density estimation algorithm to use.
#' @param Nbootstrap The number of bootstrap samples to generate.
#' @return A list containing bootstrap results and confidence intervals.
#' @examples
#' set.seed(1234); dat <- c(rnorm(200, 4), rnorm(50, 8, 2))
#' bestAlgorithm <- getBestAlgorithmByCv(dat)
#' res <- getBootstrapResult(dat, bestAlgorithm)
#' plotBootstrapResult(res)
#' 
#' @export 
getBootstrapResult <- function(dat, densFn, Nbootstrap = 400) {
  boot <- Bootstrap$new(Nbootstrap, dat, densFn, 0)
  # exec Bootstrap
  res <- boot$doAllProcess(show=FALSE)
  # get confident interval
  ci <- boot$getCI()
  return (list(bootRes = res, ci = ci, boot = boot)) 
}

#' Plot bootstrap results for density estimation.
#'
#' This function plots bootstrap results, including confidence intervals and density
#' estimates, for a density estimation algorithm.
#'
#' @param bootstrapRes The results of the bootstrap analysis.
#' @param method The plotting method (1 or 2).
#' @param showDetectedPeaks Whether to display detected peaks.
#' @examples
#' set.seed(1234); dat <- c(rnorm(200, 4), rnorm(50, 8, 2))
#' bestAlgorithm <- getBestAlgorithmByCv(dat)
#' res <- getBootstrapResult(dat, bestAlgorithm)
#' plotBootstrapResult(res, method = 1)
#' plotBootstrapResult(res, method = 2)
#' 
#' @export 
plotBootstrapResult <- function (bootstrapRes, method = 1, showDetectedPeaks = TRUE) {
  plot(
    NA,
    NA,
    xlim=range(dat),
    ylim=range(c(bootstrapRes$ci$upperCI,bootstrapRes$ci$lowerC)),
    xlab = "age",
    ylab = "density"
  )

  if (method == 1) {
    polygon(x = c(bootstrapRes$ci$x, rev(bootstrapRes$ci$x)), y = c(bootstrapRes$ci$lowerCI, rev(bootstrapRes$ci$upperCI)), col = "lightblue", border = NA)
    lines(bootstrapRes$boot$baseDens, lwd = 1.5)

  } else if (method == 2) {
    lines(bootstrapRes$boot$baseDens, lwd = 1.5)
    if (showDetectedPeaks) {
      x = bootstrapRes$boot$filteredPeaksInBootstrapResult %>% lapply(function(x) x$peaksX) %>% unlist
      y = bootstrapRes$boot$filteredPeaksInBootstrapResult %>% lapply(function(x) x$peaksY) %>% unlist
      points(x, y, pch = 20, col="grey50")
    }
    points(bootstrapRes$bootRes$mu, bootstrapRes$bootRes$y, pch = 25, col="black", cex = 2, bg = "yellowgreen")

  } else {
    stop("method is 1 or 2")
  }

}

#' Perform XBADZUPA analysis
#'
#' This function performs XBADZUPA peak detection analysis on the provided data.
#'
#' @param dat The data for peak detection.
#' @param Nbootstrap The number of bootstrap samples to generate.
#' @param ci The confidence level for peak detection.
#' @param showDetectedPeaks Whether to display detected peaks.
#' @param ignoreWarn Whether to ignore warning messages.
#' @examples
#' set.seed(1234); dat <- c(rnorm(200, 4), rnorm(50, 8, 2))
#' xbadzupa(dat, 100)
#' 
#' @export 
xbadzupa <- function(dat, Nbootstrap = 400, ci = 0.90, showDetectedPeaks = TRUE, ignoreWarn = TRUE) {

  # density estimation
  densFn <- getBestAlgorithmByCv(dat)

  # bootstrap
  if (ignoreWarn) {
    old <- options("warn")$warn
    options(warn = -1)
    res <- getBootstrapResult(dat, densFn, Nbootstrap)
    options(warn = old)
  } else {
    res <- getBootstrapResult(dat, densFn, Nbootstrap)
  }

  cat("--------XBADZUPA results of peak detection algorithm-------\n")
  cat("mu: The ages of peaks.\n")
  cat("pi: The confidence of the peaks, ranging between 0 and 1.\n")
  cat("s: The statistical error (1 sigma) in the ages.\n")
  cat("prominence *1: The prominence value of the ages.\n")
  cat("y: The probability density values of the ages.\n")
  cat("*1 Refer to the original research paper for further details.\n")
  cat("..............................................................\n")
  print(res$bootRes)
  cat("------------------------------------------------------------\n")

  split.screen(figs = c(2, 1))
  screen(1)
  plotBootstrapResult(res, 1)
  rug(dat)
  screen(2)
  plotBootstrapResult(res, 2)
  rug(dat)

}
 

