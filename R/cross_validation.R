#' R6 Class Bootstrap
#' text
Cv <- R6::R6Class(

  "cv",
  private = list(
    cvData = NULL
  ),

  public = list(

      dat = NULL,
      algorithms = NULL,
      cvN = NULL,
      watcher = NULL,
      popDens = NULL,
      progressAlgorithm = 1,

      initialize = function(algorithms, dat, cvN, watcher=function(index, all){ }) {
        self$algorithms <- algorithms
        self$dat <- dat
        self$cvN <- cvN
        self$watcher <- watcher
      },

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
      },

      rlcv = function(dat, dens, vdat) {
        a <- self$getA(dat)
        l <- sum(-self$logStar(tools$fitDistToX(sort(vdat), dens)$y, a)) / length(vdat)
        uy <- dens$y[dens$y > a] / sum(dens$y)
        ly <- (dens$y[dens$y <= a]) ^ 2 / sum(dens$y)
        b <- sum(uy) + sum(1 / (2 * a) * ly)
        if (b == 0) {
          print(paste("ly", ly))
          print(paste("uy", uy))
          stop()
        }
        print(paste('(likelihood, bias)=(', l, b, ')'), quote = FALSE)
        return(l + b)
      },

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

      #' @description
#' Change hair color.
#' @param seed New hair color.
#' @examples
#' algorithms <- list(
#'     densityDefault = function(x) density(x),
#'     denstiySJ = function(x) density(x, bw = "SJ")
#' )
#' dat <- rnorm(100)
#' cvN <- 10
#' cv <- Cv$new(algorithms, dat, cvN)
#' res <- cv$exeCVinAllAlgorithm()
      exeCVinAllAlgorithm = function(seed = 1) {
        res <- data.frame(index = c(), algorithm = c(), score = c());
        cvScores <- numeric(length(self$algorithms))
        private$cvData <- self$divideVector(self$cvN, self$dat)
        for (i in seq(self$algorithms)) {
          self$progressAlgorithm = i;
          print(paste('-------start alogrithm', names(self$algorithms)[i], "------------"), quote = FALSE)
          fn <- self$algorithms[[i]]
          cvScore <- self$exeCV(fn, self$cvN, seed, algorithmName=names(self$algorithms[i]))
          cvScores[i] <- mean(cvScore$score)
          res <- rbind(res, cvScore)
          print(paste('> mean of cv score', mean(cvScore$score)), quote = FALSE)
        }
        print("--------best algorithm----------------------------", quote = FALSE)
        print(names(self$algorithms[which.min(cvScores)]))
        print("--------------------------------------------------", quote = FALSE)
        return(res)
      }
  )
)
