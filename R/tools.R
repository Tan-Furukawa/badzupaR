Tool <- R6::R6Class(
    "tools",
    public = list(

        fitDistToX = function(x, dist) {
          tx <- dist$x
          ty <- dist$y
          y <- numeric(length(x))
          for (i in 1:length(x)) {
            j <- bsearchtools::lb(tx, x[i])
            if (j <= 1) {
              y[i] <- ty[1]
              next
            }
            if (j > length(tx)) {
              y[i] <- ty[length(ty)]
              # break
              next
            }
            y[i] <- (ty[j] - ty[j - 1]) / (tx[j] - tx[j - 1]) * (x[i] - tx[j]) + ty[j]
          }
          return(list(x = x, y = y))
        },

        sampleFromPopDist = function(n, popDist, seed = 1234) {
          set.seed(NULL);
          return(sample(popDist$x, n, prob = (popDist$y) / sum(popDist$y), replace = TRUE))
        },

        KLdistance = function(sampleDist, popDist) {
          x <- popDist$x
          modifSampleDist <- self$fitDistToX(x, sampleDist)
          q <- modifSampleDist$y
          p <- popDist$y
          dx <- diff(modifSampleDist$x)
          dx <- c(dx, dx[length(dx)])
          kl <- (popDist$y * log(p / q)) * dx
          kl[is.na(kl)] <- 0
          return(sum(kl))
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

        KLstarDistance = function(sampleDist, popDist, dat, basedOnPop = T) {
          a <- self$getA(dat)
          x <- popDist$x
          modifSampleDist <- self$fitDistToX(x, sampleDist)
          q <- modifSampleDist$y
          p <- popDist$y
          dx <- diff(modifSampleDist$x)
          dx <- c(dx, dx[length(dx)])
          p[is.na(p)] <- 0
          q[is.na(q)] <- 0
          if (basedOnPop) {
            kl <- (popDist$y * (self$logStar(p, a) - self$logStar(q, a))) * dx
          } else {
            kl <- (popDist$y * (-self$logStar(q, a))) * dx
          }
          kl[is.na(kl)] <- 0
          return(sum(kl))
        },

        L2distance = function(sampleDist, popDist) {
          x <- popDist$x
          modifSampleDist <- self$fitDistToX(x, sampleDist)
          dx <- diff(modifSampleDist$x)
          dx <- c(dx, dx[length(dx)])
          l2 <- sum((modifSampleDist$y - popDist$y) ^ 2 * dx)
          return(l2)
        },

        addRandomNoise = function(vec, noise) {
          newVec <- numeric(length(vec))
          for (i in seq(vec)) {
            newVec[i] <- vec[i] + rnorm(1, 0, noise)
          }
          return(newVec)
        },

        getTime = function() {
          return(as.numeric(as.POSIXct(Sys.time())))
        },

        normalizeDat = function(dat) {
          return(dat - min(dat)) / (max(dat) - min(dat))
        },

        UnNormalizeDat = function(ndat, dat) {
          return(ndat * (max(dat) - min(dat)) + min(dat))
        },

        estimateCoodinateFromThreePoint = function(x, y) {
          if (length(x) != 3) stop("x must 3 element")
          if (length(y) != 3) stop("x must 3 element")
          if ((y[1] <= y[2] & y[2] > y[3]) | (y[1] < y[2] & y[2] >= y[3])) {
            d = x[3] - x[2];
            p = (y[3] + y[1] - 2.0 * y[2]) / (2.0 * d * d);
            q = (y[3] - y[1]) / (2.0 * d) - 2.0 * p * x[2];
            r = y[2] - p * x[2] * x[2] - q * x[2];
            xx = -q / (2.0 * p);
            yy = p * xx * xx + q * xx + r;
            return(list(x = xx, y = yy))
          } else {
            return(NULL)
          }
        },

        applyWithoutDrop = function(X, MARGIN, FUN, ...) {
          if (1 %in% dim(X)) {
            return(t(apply(t(X), 2, FUN, ...)))
          } else {
            return(apply(X, MARGIN, FUN, ...))
          }
        },

        makeUniformNoiseFromData = function(dat) {
          sdat <- sort(dat)
          ddat <- diff(sort(sdat))
          u <- sdat + c(ddat / 2, 0)
          l <- sdat - c(0, ddat / 2)
          res <- sample(1:length(dat), replace = T) %>% 
            sapply(function(i) runif(1, min=l[i], max=u[i]))
          return(res)
        }

        )
)

tools <- Tool$new()
