#' Cross Validation for xbadzupaR
#'
#' This R6 class provides cross-validation functionality for evaluating algorithms.
#' 
#' @example 
#' algorithms <- list(
#'     denstiyBOTEV = function(x) IsoplotR::kde(x, plot=FALSE),
#'     densityADEBA = density.adeba
#' )
#' dat <- rnorm(100)
#' cv <- Cv$new(algorithms, dat)
#' cv$exeCVinAllAlgorithm()
#'
#' @importFrom R6 R6Class
#' @export
Cv <- R6::R6Class(

  "cv",
  private = list(
    cvData = NULL,
    divideVector = function(N, vec) {
      vecRem <- 1:length(vec)
      distribute <- length(vec) %% N
      eachLen <- length(vec) %/% N
      validation <- list()

      for (i in 1:N) {
        if (i <= distribute) {
          to <- eachLen + 1
        } else {
          to <- eachLen
        }
        validation[[i]] <- vecRem[1:to]
        vecRem <- vecRem[(to + 1):length(vecRem)]
      }

      res <- list()
      for (i in seq(validation)) {
        res[[i]] <- list(
              train = vec[-validation[[i]]],
              validation = vec[validation[[i]]]
          )
      }
      return(res)
    },

    logStar = function(x, a) {
      res <- numeric(length(x))
      for (i in 1:length(x)) {
        if (x[i] > a) {
          res[i] <- (log(x[i]))
        } else {
          res[i] <- (log(a) - 1 + x[i] / a)
        }
      }
      return(res)
    },

    getA = function(dat) {
      n <- length(dat)
      return(sd(dat) ^ (-1) * (2 * pi) ^ (-0.5) * gamma(0.5) * (log(n)) ^ (0.5) / n)
    }
  ),

  public = list(

      cvN = 10,
      dat = NULL,
      algorithms = NULL,
      watcher = NULL,
      popDens = NULL,
      usingLogStar = TRUE,
      progressAlgorithm = 1,

      #' Initialize the Cv class
      #'
      #' @param algorithms A list of algorithms to be evaluated.
      #' @param dat The input data for cross-validation.
      #' @param watcher A callback function to monitor progress.
      #'
      #' @return A Cv object.
      initialize = function(algorithms, dat, watcher=function(index, all){ }) {
        self$algorithms <- algorithms
        self$dat <- dat
        self$cvN <- self$cvN
        self$watcher <- watcher
      },

      #' Calculate Log-Likelihood with Bias for Cross-Validation
      #'
      #' @param dat A numeric vector representing the input data used for density estimation.
      #' @param dens A "density" object representing the density estimation.
      #' @param vdat A numeric vector representing the validation data.
      #' @return The log-likelihood with added bias.
      #'
      rlcv = function(dat, dens, vdat) {
        if (self$usingLogStar) {
          a <- private$getA(dat)
          l <- sum(-private$logStar(tools$fitDistToX(sort(vdat), dens)$y, a)) / length(vdat)
          uy <- dens$y[dens$y > a] / sum(dens$y)
          ly <- (dens$y[dens$y <= a])^2 / sum(dens$y)
          b <- sum(uy) + sum(1 / (2 * a) * ly)
          
          if (b == 0) {
            cat(paste("ly", ly, "\n"))
            cat(paste("uy", uy, "\n"))
            stop()
          }
          cat(paste('(likelihood, bias)=(', l, b, ')\n'))
          return(l + b)
        } else {
          l <- sum(-log(tools$fitDistToX(sort(vdat), dens)$y)) / length(vdat)
          cat(paste('(likelihood, bias)=(', l, ', NA)\n'))
          return(l)
        }

      },

      #' Execute cross-validation for a specific algorithm
      #'
      #' @param fn The algorithm function to be evaluated.
      #' @param seed A seed value for reproducibility.
      #' @param save A logical indicating whether to save results.
      #' @param algorithmName The name of the algorithm.
      #'
      #' @return A data frame containing cross-validation results.
      exeCV = function(fn, seed = NULL, save = TRUE, algorithmName="") {
        scores <- c()
        res = data.frame(index = c(), algorithm = c(), score = c())
        for (i in 1:self$cvN) {
          self$watcher (i + self$cvN * (self$progressAlgorithm - 1), self$cvN * length(self$algorithms))
          vdat <- private$cvData[[i]]$validation

          dens <- fn(private$cvData[[i]]$train)
          score <- self$rlcv(self$dat, dens, vdat)
          res <- rbind(res, data.frame(index = i, algorithm = algorithmName, score = score))
        }
        return(res)
      },

      #' Execute cross-validation for all algorithms
      #'
      #' @param seed A seed value for reproducibility.
      #' 
      #' @return A data frame containing cross-validation results for all algorithms.
      #' 
      #' @examples
      #'  algorithms <- list(
      #'      denstiyBOTEV = function(x) IsoplotR::kde(x, plot=FALSE),
      #'      densityADEBA = density.adeba
      #'  )
      #'  dat <- rnorm(100)

      #'  cv <- Cv$new(algorithms, dat)
      #'  cv$exeCVinAllAlgorithm()
      #'  
      exeCVinAllAlgorithm = function(seed = 1) {
        res <- data.frame(index = c(), algorithm = c(), score = c());
        cvScores <- numeric(length(self$algorithms))
        private$cvData <- private$divideVector(self$cvN, self$dat)
        for (i in seq(self$algorithms)) {
          self$progressAlgorithm = i;
          cat(paste('-------start alogrithm', names(self$algorithms)[i], "------------\n"))
          fn <- self$algorithms[[i]]
          cvScore <- self$exeCV(fn, self$cvN, seed, algorithmName=names(self$algorithms[i]))
          cvScores[i] <- mean(cvScore$score)
          res <- rbind(res, cvScore)
          cat(paste('mean of cv score', mean(cvScore$score),"\n"))
        }
        cat("--------best algorithm----------------------------\n")
        cat(paste(names(self$algorithms[which.min(cvScores)])), "\n")
        cat("--------------------------------------------------\n")
        return(res)
      }
  )
)
