# 正規分布＋一様ノイズでフィッティング
# ====================================

makeLikelihood <- function(mu0, c) {
  return(function(pi1, pi2, sd) {
    return(function(x) {
        res <- matrix(
              data = c(
                c(pi1 * dnorm(x, mean = mu0, sd = sd)),
                c(pi2 * c + numeric(length(x)))
              ), ncol=2, byrow=F
          )
        colnames(res) <- c('norm','noise')
        return (res)
    })
  })
}
# print(makeLikelihood(0, 1)(.5,.5, 1)(1:10*0.1))

getSd = function (dat, likelihood, mu0) {
    resSum <- sum( likelihood[,1] )
    resSd <- sqrt( sum( likelihood[,1] * (dat - mu0)^2 ) / resSum )
    return (resSd)
}
# print(getSd((1:10)*0.1, makeLikelihood(0, 1)(.5,.5, 1)((1:10)*0.1), 0))

applyWithoutDrop = function(X, MARGIN, FUN, ...) {
  if (1 %in% dim(X)) {
    return(t(apply(t(X), 2, FUN, ...)))
  } else {
    return(apply(X, MARGIN, FUN, ...))
  }
}

normalizeLikelihood = function (likelihood) {
    I = ncol(likelihood)
    res <- t(applyWithoutDrop(likelihood, 1, function (x) {
        return (x / sum(x))
    }))
    return (res)
}
# print(normalizeLikelihood(makeLikelihood(0, 1)(.5,.5, 1)((1:10)*0.1)))

getPi = function (nLikelihood) {
    if (is.null(dim(nLikelihood))) {
        return (1)
    } else {
        s <- colSums(nLikelihood)
        return (s / sum(s))
    }
}
# res <- normalizeLikelihood(makeLikelihood(0, 1)(.5,.5, 1)((1:10)*0.1))
# print(getPi(res))
        
doEM = function (dat, sd0, mu0, pi1_0, pi2_0, showProgress=TRUE) {
    pi1 <- pi1_0
    pi2 <- pi2_0
    sd <- sd0
    c <- 0.1 / (max(dat) - min(dat))
    for (j in 1:100) {
        nlikelihoodList <- makeLikelihood(mu0,c)(pi1, pi2,sd)(dat) %>% normalizeLikelihood()
        sd <- getSd (dat, nlikelihoodList, mu0)
        piRes <- getPi (nlikelihoodList)
        pi1 <- piRes[1]
        pi2 <- piRes[2]
        if (showProgress) {
          # print(sd)
          # print(piRes)
            cat("GMM Process["); cat(j); cat(";sd:"); cat(sd); cat(";pi:"); cat(piRes); cat("\n")
        }
    }
    correctRatio = sum((nlikelihoodList >= 0.5)[,1]) / nrow(nlikelihoodList)
    return (list(mu=mu0, s=sd, pi=correctRatio, pi1 = pi1, pi2=pi2))
}

# ##1
# x <- seq(0,1000,length=100)
# dat <- c(rnorm(100,200,100), runif(100, min=0, max=2000))
# correct <- 0.5 * dnorm(x, 200, 100) + 0.5 * dunif(x,min=0,max=1000)
# estimate <- function(res) res$pi1 * dnorm(x, res$mu, res$sd) + res$pi2 * 1 / (max(dat) - min(dat))

# res <- doEM(dat, 10, 200, 0.5, 0.5)
# plot( x, correct, type='l',)
# rug(dat)
# lines( x,estimate(res), col='red')
# print(res)

##2
# x <- seq(0,1000,length=100)
# dat <- c(rnorm(400,200,100))
# correct <- 1 * dnorm(x, 200, 100)
# estimate <- function(res) res$pi1 * dnorm(x, res$mu, res$sd) + res$pi2 * 1 / (max(dat) - min(dat))
# res <- doEM(dat, 100, 200, 0.5, 0.5)
# print(res)
# plot( x, correct, type='l')
# rug(dat)
# lines( x,estimate(res), col='red')

## 3
# x <- seq(0,2000,length=100)
# dat <- c(rnorm(200,200,100), rnorm(200,1000,100))
# correct <- 1 * 0.5 * dnorm(x, 200, 100) + 0.5 * dnorm(x, 1000, 100)
# estimate <- function(res) res$pi1 * dnorm(x, res$mu, res$sd) + res$pi2 * 1 / (max(dat) - min(dat))
# res <- doEM(dat, 200, 200, 0.5, 0.5)
# print(res)
# plot( x, correct, type='l')
# rug(dat)
# lines( x,estimate(res), col='red')
