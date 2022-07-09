Gmm <- R6::R6Class(
    "gmm",
    public = list(
        makeLikelihoodList = function (dat, s, mu, pi) {
            I <- length(s)
            n <- length(dat)
            res <- matrix(NA, nrow=n, ncol = I)
            for (j in 1:n) {
                for (i in 1:I) {
                    res[j,i] <- pi[i] * dnorm(dat[j], mu[i], s[i])
                }
            }
            return(res)
        },
        
        getMu = function (dat, likelihoodList, mu) {
            # fixMu <- paramList$fixMu
            # 
            # I <- ncol(likelihoodList)
            # j <- which.min(abs(mu0 - fixMu))
            # mu <- numeric(length(mu0))
            # sum <- numeric(length(mu0))
            # for (i in 1:I) {
            #   if (i == j) {
            #     mu[j] = fixMu
            #   } else {
            #     sum[i] <- sum(likelihoodList[,i])
            #     mu[i] <- sum(likelihoodList[,i] * data)
            #   }
            # }
            return (mu)
        },
        
        getS = function (dat, likelihoodList, mu) {
            I = ncol(likelihoodList)
            resS <- numeric(I)
            resSum <- numeric(I)
            for (i in 1:I) {
                resSum[i] <- sum( likelihoodList[,i] )
                resS[i] <- sqrt( sum( likelihoodList[,i] * (dat - mu[i])^2 ) / resSum[i] )
            }
            return (resS)
        },
        
        normalizeLikelihoodList = function (likelihoodList) {
            I = ncol(likelihoodList)
            res <- t(tools$applyWithoutDrop(likelihoodList, 1, function (x) {
                return (x / sum(x))
            }))
            # v <- apply(likelihoodList, 1, function(x) sum(x))
            # res <- likelihoodList / t(matrix(v, nrow=ncol(likelihoodList), ncol=length(v), byrow=TRUE))
            return (res)
        },
        
        getPi = function (nlikelihoodList) {
            if (is.null(dim(nlikelihoodList))) {
                return (1)
            } else {
                s <- colSums(nlikelihoodList)
                return (s / sum(s))
            }
        },
        
        doGMM = function (dat, s, mu, pi, showProgress=TRUE) {
            for (j in 1:10) {
                nlikelihoodList <- self$makeLikelihoodList(dat,s,mu,pi) %>% self$normalizeLikelihoodList()
                mu <- self$getMu (dat, nlikelihoodList, mu)
                s <- self$getS (dat, nlikelihoodList, mu)
                pi <- self$getPi (nlikelihoodList)
                if (showProgress) {
                    cat("GMM Process["); cat(j);cat("]=> mu:"); cat(mu); cat(";sd:"); cat(s); cat(";pi:"); cat(pi); cat("\n")
                }
            }
            return (list(mu=mu, s=s, pi=pi))
        }
    )
)

# dat <- rnorm(50)
gmm <- Gmm$new()
# dat <- c(dat,rnorm(100,3))
# s0 <- c(1)
# mu0 <- c(0)
# pi0 <- c(0.5)
# x <- seq(-10,10,length=100)
# res<- gmm$doGMM(dat,s0,mu0,pi0)
# plot(density(dat),type="l")
# lines(x,
#       res$pi[1] * dnorm(x,res$mu[1],res$s[1])
#       )
# rug(dat)

# gmm <- gmm$new()
# s0 <- ini$iniSd
# mu0 <- ini$iniMu
# pi0 <- ini$iniSd
# x <- seq(0,3000,length=100)
# res <- gmm$doGMM(dd,s0,mu0,pi0)
# plot(density(dd),type="l")
# lines(x,dnorm(x,res$mu,res$s))
# rug(dd)

