#' R6 Class Bootstrap
#'
#' This R6 class, Bootstrap, is designed for bootstrapping density estimation and peak detection.
#' It allows you to perform bootstrap resampling, detect peaks in the bootstrapped results, and
#' estimate Gaussian Mixture Models (GMM) for the detected peaks.
#'
#' @usage
#' set.seed(123); dat <- rnorm(100)
#' densFn <- density
#' Nbootstrap <- 100
#' noiseSize <- 0
#' boot <- Bootstrap$new(Nbootstrap, dat, densFn, noiseSize)
#' # exec Bootstrap
#' res <- boot$doAllProcess(show=FALSE)
#' # get confident interval
#' ci <- boot$getCI()
#' # plot result
#' plot(NA,NA,xlim=c(-3,3),ylim=c(0,0.5))
#' polygon(x = c(ci$x, rev(ci$x)), y = c(ci$lowerCI, rev(ci$upperCI)), col = "lightblue", border = NA)
#' lines(boot$baseDens)
#' curve(dnorm,-3,3,add=TRUE,col="red")
#'
#' @section Public Fields:
#' \describe{
#'   \item{watcher}{A function to monitor the progress of bootstrapping.}
#'   \item{Nbootstrap}{The number of bootstrap samples to generate.}
#'   \item{dat}{The original data for density estimation.}
#'   \item{ci}{A numeric vector representing confidence intervals.}
#'   \item{noiseSize}{The standard deviation of random noise to add during bootstrapping.}
#'   \item{densFn}{A function for estimating the density of the data.}
#'   \item{baseDens}{The density estimate of the original data.}
#'   \item{bootstrapResult}{A matrix containing bootstrap results.}
#'   \item{peaksInBootstrapResult}{A list of peak information from the bootstrap results.}
#'   \item{filteredPeaksInBootstrapResult}{Filtered peak information after applying thresholds.}
#'   \item{prominenceThreshold}{Threshold for peak prominence.}
#'   \item{relativeHeight}{Threshold for relative peak height.}
#'   \item{bootstrapAllPeaksX}{A matrix containing clustered peaks from the base density.}
#'   \item{basePeaks}{Filtered base peaks.}
#'   \item{method}{The method used for bootstrapping ('naive' or 'smooth').}
#' }
#'
#' @section Public Methods:
#' \describe{
#'   \item{initialize}{Initialize the Bootstrap object with the specified parameters.}
#'   \item{getBaseDens}{Estimate the density of the original data.}
#'   \item{exeBootstrap}{Perform bootstrap resampling.}
#'   \item{plotAllBootResult}{Plot all bootstrapped results.}
#'   \item{detectAllPeaks}{Detect peaks in a density estimate.}
#'   \item{detectPeaks}{Detect peaks based on specified criteria.}
#'   \item{detectPeaksInBootstrapResult}{Detect peaks in the bootstrap results.}
#'   \item{filterProminence}{Filter peaks based on prominence and relative height thresholds.}
#'   \item{filterAllProminence}{Filter all peaks in the bootstrap results.}
#'   \item{clusterPeaks}{Cluster peaks around a target X value.}
#'   \item{clusterAllPeaks}{Cluster all base peaks.}
#'   \item{checkBootstrapAllPeaks}{Plot density of clustered peaks.}
#'   \item{makeInitialGmmParam}{Generate initial GMM parameters for a specific bootstrap result.}
#'   \item{execGmm}{Perform Gaussian Mixture Model (GMM) estimation for a specific bootstrap result.}
#'   \item{execAllGmm}{Perform GMM estimation for all clustered peaks.}
#'   \item{getCI}{Get confidence intervals for density estimation.}
#'   \item{doAllProcess}{Execute the entire Bootstrap process.}
#' }
#'
#' @importFrom R6 R6Class
#' @export
Bootstrap <- R6::R6Class(

    "bootstrap",

    public = list(
        watcher = NULL,
        Nbootstrap = NA,
        dat = NA,
        ci = NA,
        noiseSize = NA,
        densFn = NULL,
        baseDens = NULL,
        bootstrapResult = NULL,
        peaksInBootstrapResult = NULL,
        filteredPeaksInBootstrapResult = NULL,
        prominenceThreshold = 0.01,
        relativeHeight = 0.00,
        # prominenceThresholdAtBootstrapPeaksDetection = 0.05,
        # relativeHeightAtBootstrapPeaksDetection = 0.05,
        bootstrapAllPeaksX = NULL,
        basePeaks = NULL,
        basePeaksMean = NULL,
        # method = 'empirical',
        method = 'smooth',

        #' initialize
        #'
        #' This function initializes the Bootstrap object with essential parameters and functions
        #' for bootstrapping and density estimation.
        #'
        #' @param Nbootstrap The number of bootstrap samples to generate.
        #' @param dat The original dataset for bootstrapping.
        #' @param densFn The density estimation function for the dataset.
        #' @param noiseSize The amount of random noise to add to resampled data.
        #' @param watcherFn A custom function to monitor the progress of bootstrapping (default is a no-op function).
        initialize = function(Nbootstrap, dat, densFn, noiseSize, watcherFn = (function(i, N) { })) {
          self$Nbootstrap <- Nbootstrap
          self$dat <- dat
          self$densFn <- densFn
          self$noiseSize <- noiseSize
          self$watcher <- watcherFn
        },

        #' getBaseDens
        #'
        #' This function computes the density function of the original dataset using the provided density estimation function.
        #'
        getBaseDens = function() {
          self$baseDens <- self$densFn(self$dat)
        },

        #' exeBootstrap
        #'
        #' This function performs bootstrapping to generate multiple resampled datasets
        #' and estimates their density functions. It is used to create the bootstrap result
        #' for further analysis.
        #'
        #' @param progress A logical value indicating whether to display progress during bootstrapping (default is TRUE).
        #'
        exeBootstrap = function(progress = TRUE) {
          m <- length(self$baseDens$x)
          pb <- NULL
          if (progress) {
            pb <- progress::progress_bar$new(format = "progress of bootstrap [:bar] :current/:total (:percent)", total = self$Nbootstrap)
          }
          self$bootstrapResult <- matrix(NA, nrow = self$Nbootstrap, ncol = m)

          for (i in 1:self$Nbootstrap) {
            if (progress) pb$tick()
            self$watcher(i, self$Nbootstrap)
            if (self$method == 'smooth') {
              sampledDat <- tools$sampleFromPopDist(length(self$dat), self$baseDens, seed = i)
              sampledDat <- tools$addRandomNoise(sampledDat, self$noiseSize)
            } else if (self$method == 'empirical') {
              sampledDat <- tools$makeUniformNoiseFromData(dat)
            }
            else {
              sampledDat <- sample(self$dat, length(self$dat), replace = TRUE)
              sampledDat <- tools$addRandomNoise(sampledDat, self$noiseSize)
            }

            sampledDens <- self$densFn(sampledDat)
            self$bootstrapResult[i,] <- tools$fitDistToX(self$baseDens$x, sampledDens)$y
          }
        },

        #' plotAllBootResult
        #'
        #' This function generates a plot to visualize the bootstrap results and the original dataset's density function.
        # It overlays lines representing each bootstrapped density function with transparency to show the variation.
        #'
        plotAllBootResult = function() {
          for (i in 1:self$Nbootstrap) {
            lines(self$baseDens$x, self$bootstrapResult[i,], col = rgb(0, 0, 1, 0.1))
          }
          lines(self$baseDens)
        },

        detectAllPeaks = function(resDens) {
          y <- resDens$y
          x <- resDens$x
          peaks <- c()
          for (i in 2:(length(x) - 1)) {
            if (y[i - 1] < y[i] & y[i] > y[i + 1]) {
              peaks <- c(peaks, i)
            }
          }
          return(peaks)
        },

        detectPeaks = function(resDens) {
          peaks <- self$detectAllPeaks(resDens)
          prominence <- numeric(length(peaks))
          width <- numeric(length(peaks))

          y <- resDens$y
          x <- resDens$x
          n <- length(y)

          if (length(peaks) == 0) {
            warning("no peaks")
            return(numeric(0))
          }

          maxPeaksY <- max(y[peaks])

          for (j in 1:length(peaks)) {
            i <- peaks[j]
            maxI <- length(y) - 1
            minI <- 1
            for (k in i:n) {
              if (y[i] < y[k]) {
                maxI <- k
                break
              }
            }

            for (k in i:1) {
              if (y[i] < y[k]) {
                minI <- k
                break
              }
            }

            largerAreaMinY <- min(y[i:maxI])
            smallerAreaMinY <- min(y[minI:i])

            if (largerAreaMinY > smallerAreaMinY) {
              prominence[j] <- y[i] - largerAreaMinY
            } else {
              prominence[j] <- y[i] - smallerAreaMinY
            }
          }

          peaksX = x[peaks]
          peaksY = y[peaks]
          for (i in 1:length(peaks)) {
            j = peaks[i]
            if (j == 1) {
              #do nothing
            } else if (j == length(y)) {
              #do nothing
            } else {
              modifXY = tools$estimateCoodinateFromThreePoint(x = c(x[j - 1], x[j], x[j + 1]), y = c(y[j - 1], y[j], y[j + 1]))
              peaksX[i] = modifXY$x
              peaksY[i] = modifXY$y
            }
          }

          return(list(
                prominence = prominence,
                rprominence = prominence / y[peaks],
                peaksY = peaksY,
                peaksX = peaksX
            ))
        },

        detectPeaksInBootstrapResult = function() {
          res <- vector(mode = "list", length = self$Nbootstrap)
          for (i in 1:self$Nbootstrap) {
            res[[i]] <- self$detectPeaks(
                    list(
                        x = self$baseDens$x,
                        y = self$bootstrapResult[i,]
                    )
                )
          }
          self$peaksInBootstrapResult <- res
        },

        filterProminence = function(eachPeak) {
          useIt <- (eachPeak$peaksY > max(eachPeak$peaksY) * self$relativeHeight) & (eachPeak$rprominence > self$prominenceThreshold)
          return(list(
                prominence = eachPeak$prominence[useIt],
                rprominence = eachPeak$rprominence[useIt],
                peaksY = eachPeak$peaksY[useIt],
                peaksX = eachPeak$peaksX[useIt]
            ))
        },

        filterAllProminence = function() {
          self$filteredPeaksInBootstrapResult <- self$peaksInBootstrapResult
          for (i in 1:length(self$peaksInBootstrapResult)) {
            eachPeak <- self$peaksInBootstrapResult[[i]]
            self$filteredPeaksInBootstrapResult[[i]] <-
                    self$filterProminence(eachPeak)
          }
        },

        clusterPeaks = function(targetX) {
          eachPeakCluster <- sapply(self$filteredPeaksInBootstrapResult, function(l) {
            return(l$peaksX[which.min(abs(l$peaksX - targetX))])
          })
          return(eachPeakCluster)
        },

        clusterAllPeaks = function() {
          notFilteredBasePeaks <- self$detectPeaks(self$baseDens)
          self$basePeaks <- self$filterProminence(notFilteredBasePeaks)
          self$bootstrapAllPeaksX <- apply(as.matrix(self$basePeaks$peaksX), 1, function(x) {
            return(self$clusterPeaks(x))
          })
        },

        checkBootstrapAllPeaks = function(i) {
          plot(density(self$bootstrapAllPeaksX[, i]))
          rug(self$bootstrapAllPeaksX[, i])
        },

        makeInitialGmmParam = function(d, i) {

          basePeaksX = self$basePeaks$peaksX[i]

          dens <- density(d, bw = "sj-dpi")
          peaks <- self$detectPeaks(dens)

          escProminenceThreshold <- self$prominenceThreshold
          escRelattiveHeight <- self$relativeHeight

          # self$prominenceThreshold <- self$prominenceThresholdAtBootstrapPeaksDetection
          # self$relativeHeight <- self$relativeHeightAtBootstrapPeaksDetection

          peaks <- self$filterProminence(peaks)

          self$prominenceThreshold <- escProminenceThreshold
          self$relativeHeight <- escRelattiveHeight

          group <- numeric(length(d))
          iniSd <- numeric(length(peaks$peaksX))

          for (j in 1:length(d)) {
            group[j] <- which.min(abs(peaks$peaksX - d[j]))
          }
          for (j in 1:length(peaks$peaksX)) {
            iniSd[j] <- sd(d[group == j])
          }
          # validation 
          useIt <- !is.na(iniSd)
          iniMu <- (peaks$peaksX)[useIt]
          isBase <- which.min(abs(basePeaksX - iniMu))
          iniMu[isBase] <- basePeaksX
          return(
                list(
                    iniMu = (peaks$peaksX)[useIt],
                    iniPi = (peaks$prominence / sum(peaks$prominence))[useIt],
                    iniSd = iniSd[useIt],
                    isBase = isBase
                )
            )
        },

        execGmm = function(i) {
          d <- self$bootstrapAllPeaksX[, i]
          initialParams <- self$makeInitialGmmParam(d, i)
          isBase = initialParams$isBase
          # print(initialParams)
        #   res <- gmm$doGMM(d, initialParams$iniSd, initialParams$iniMu, initialParams$iniPi, showProgress = FALSE)
          res <- doEM(
              d, 
              initialParams$iniSd[isBase],
              initialParams$iniMu[isBase],
              0.7,
              0.3,
              showProgress = FALSE
              )
          return(res)
        },

        #' Estimate Gaussian Mixture Models (GMM) for all clustered peaks.
        #'
        #' This function iteratively estimates GMM parameters (mu, pi, and s) for each
        #' clustered peak using the result of the `execGmm` function.
        #'
        #' @return A data frame with GMM parameters for each clustered peak, including
        #' mu (mean), pi (mixture weight), s (standard deviation), prominence, and y values.
        #'
        execAllGmm = function() {
          peaksX = self$basePeaks$peaksX
          prominence = self$basePeaks$prominence
          I = length(peaksX)
          pi = numeric(I)
          s = numeric(I)
          mu = numeric(I)
          for (i in 1:I) {
            resGmm <- self$execGmm(i)
            # whichMuNearest <- which.min(abs(resGmm$mu - peaksX[i]))
            # mu[i] = resGmm$mu[whichMuNearest]
            # pi[i] = resGmm$pi[whichMuNearest]
            # s[i] = resGmm$s[whichMuNearest]
            mu[i] = resGmm$mu
            pi[i] = resGmm$pi
            s[i] = resGmm$s
          }
          return(data.frame(mu = mu, pi = pi, s = s, prominence = prominence, y = self$basePeaks$peaksY))
        },

        #' @description
        #' get confident interval
        #' @param lowerP confident level
        #' @param upper confident level
        #' @examples
        #' set.seed(123); dat <- rnorm(100)
        #' densFn <- density
        #' Nbootstrap <- 100
        #' noiseSize <- 0
        #' boot <- Bootstrap$new(Nbootstrap, dat, densFn, noiseSize)
        #' # exec Bootstrap
        #' res <- boot$doAllProcess(show=FALSE)
        #' # get confident interval
        #' ci <- boot$getCI()
        #' # plot result
        #' plot(NA,NA,xlim=c(-3,3),ylim=c(0,0.5))
        #' polygon(x = c(ci$x, rev(ci$x)), y = c(ci$lowerCI, rev(ci$upperCI)), col = "lightblue", border = NA)
        #' lines(boot$baseDens)
        #' curve(dnorm,-3,3,add=TRUE,col="red")
        #'
        getCI = function(lowerP = 0.05, upperP = 0.95) {
          res <- list(
                x = self$baseDens$x,
                lowerCI = apply(self$bootstrapResult, 2, function(x) quantile(x, probs = upperP)),
                upperCI = apply(self$bootstrapResult, 2, function(x) quantile(x, probs = lowerP))
            )
          return(res)
        },

        #' @description
        #' Executes the entire Bootstrap process, including density estimation, bootstrapping, peak detection,
        #' and Gaussian Mixture Model (GMM) estimation for clustered peaks.
        #'
        #' @param show A logical value indicating whether to display the results (default is TRUE).
        #' @param progress A logical value indicating whether to display progress bars during bootstrapping (default is TRUE).
        #'
        #' @return A data frame containing GMM parameters for the clustered peaks.
        #'
        #' @examples
        #' set.seed(123); dat <- rnorm(100)
        #' densFn <- density
        #' Nbootstrap <- 100
        #' noiseSize <- 0
        #' boot <- Bootstrap$new(Nbootstrap, dat, densFn, noiseSize)
        #' # exec Bootstrap
        #' res <- boot$doAllProcess(show=FALSE)
        #' print(res)
        #'
        doAllProcess = function(show = TRUE, progress = TRUE) {
          self$getBaseDens()
          self$exeBootstrap(progress)
          if (show) {
            self$plotAllBootResult()
          }
          self$detectPeaksInBootstrapResult()
          self$filterAllProminence()
          self$clusterAllPeaks()
          res <- self$execAllGmm()
          # ci <- self$getCI()
          return(res)
        }
    )
)
