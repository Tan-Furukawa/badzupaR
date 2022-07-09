#' R6 Class Bootstrap
#' text
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
        relativeHeight = 0.01,
        prominenceThresholdAtBootstrapPeaksDetection = 0.05,
        relativeHeightAtBootstrapPeaksDetection = 0.05,
        bootstrapAllPeaksX = NULL,
        basePeaks = NULL,
        method = 'naive',
        
        initialize = function(Nbootstrap, dat, densFn, noiseSize, watcherFn = (function(i,N){})) {
            self$Nbootstrap <- Nbootstrap
            self$dat <- dat
            self$densFn <- densFn
            self$noiseSize <- noiseSize
            self$watcher <- watcherFn
        },
        
        getBaseDens = function() {
            self$baseDens <- self$densFn(self$dat)
        },
        
        exeBootstrap = function(progress=TRUE) {
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
                } else {
                    sampledDat <- sample(self$dat, length(self$dat), replace=TRUE)
                    sampledDat <- tools$addRandomNoise(sampledDat, self$noiseSize)
                }
               
                # sampledDat <- sample(sdat, length(sdat), replace=T)
                # sampledDat <- addRandomNoise(sampledDat, 2)
                sampledDens <- self$densFn(sampledDat)
                self$bootstrapResult[i, ] <- tools$fitDistToX(self$baseDens$x, sampledDens)$y
            }
        },
        
        plotAllBootResult = function() {
            # plot(self$baseDens, type = "l")
            for (i in 1:self$Nbootstrap) {
                lines(self$baseDens$x, self$bootstrapResult[i, ], col = rgb(0, 0, 1, 0.1))
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
                } else if (j == length(y)){
                    #do nothing
                } else {
                    modifXY = tools$estimateCoodinateFromThreePoint(x=c(x[j-1],x[j],x[j+1]),y=c(y[j-1],y[j],y[j+1]))
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
                        y = self$bootstrapResult[i, ]
                    )
                )
            }
            self$peaksInBootstrapResult <- res
        },
        
        filterProminence = function(eachPeak) {
            useIt <-  (eachPeak$peaksY > max(eachPeak$peaksY) * self$relativeHeight) & (eachPeak$rprominence > self$prominenceThreshold)
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
                return (self$clusterPeaks(x))
            })
        },
        
        checkBootstrapAllPeaks = function (i) {
            plot(density(self$bootstrapAllPeaksX[,i]))
            rug(self$bootstrapAllPeaksX[,i])
        },
        
        makeInitialGmmParam = function (d, i) {
            
            basePeaksX = self$basePeaks$peaksX[i]
            
            dens <- density(d, bw="sj-dpi")
            peaks <- self$detectPeaks(dens)
            
            escProminenceThreshold <- self$prominenceThreshold
            escRelattiveHeight <- self$relativeHeight
            
            self$prominenceThreshold <-  self$prominenceThresholdAtBootstrapPeaksDetection
            self$relativeHeight <- self$relativeHeightAtBootstrapPeaksDetection
            
            peaks <- self$filterProminence(peaks)
            
            self$prominenceThreshold <- escProminenceThreshold
            self$relativeHeight <- escRelattiveHeight
            
            
            group <- numeric(length(d))
            iniSd <- numeric(length(peaks$peaksX))
            
            for (j in 1:length(d)) {
                group[j] <- which.min(abs(peaks$peaksX - d[j]))
            }
            for (j in 1:length(peaks$peaksX)) {
                iniSd[j] <- sd (d[group == j] )
            }
            # validation 
            useIt <- !is.na(iniSd)
            iniMu <- (peaks$peaksX)[useIt]
            isBase <- iniMu[which.min(abs(basePeaksX - iniMu))]
            iniMu[isBase] <- basePeaksX
            return (
                list(
                    iniMu = (peaks$peaksX)[useIt],
                    iniPi =  (peaks$prominence/ sum(peaks$prominence))[useIt],
                    iniSd = iniSd[useIt]
                )
            )
        },
        
        execGmm = function (i) {
            d <- self$bootstrapAllPeaksX[,i]
            initialParams <- self$makeInitialGmmParam(d,i)
            res <- gmm$doGMM(d, initialParams$iniSd, initialParams$iniMu, initialParams$iniPi, showProgress=FALSE)
            return (res)
        },
        
        execAllGmm = function () {
            peaksX = self$basePeaks$peaksX
            prominence = self$basePeaks$prominence
            I = length(peaksX)
            pi = numeric(I)
            s = numeric(I)
            mu = numeric(I)
            for (i in 1:I) {
                resGmm <- self$execGmm(i)
                whichMuNearest <- which.min(abs(resGmm$mu - peaksX[i]))
                mu[i] = resGmm$mu[whichMuNearest]
                pi[i] = resGmm$pi[whichMuNearest]
                s[i] = resGmm$s[whichMuNearest]
            }
            return (data.frame(mu = mu, pi = pi, s = s, prominence=prominence, y=self$basePeaks$peaksY))
        },
        #' @description
        #' Change hair color.
        #' @param show New hair color.
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
        getCI = function (lowerP = 0.05, upperP = 0.95) {
            res <- list(
                x = self$baseDens$x,
                lowerCI = apply(self$bootstrapResult, 2, function (x) quantile(x, probs=upperP)),
                upperCI = apply(self$bootstrapResult, 2, function (x) quantile(x, probs=lowerP))
            )
            return(res)
        },
        
        #' @description
        #' Change hair color.
        #' @param show New hair color.
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
        doAllProcess = function (show=TRUE, progress=TRUE) {
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
            return (res)
        }
    )
)

