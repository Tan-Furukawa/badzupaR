library(magrittr)

densityFn <- function(x) IsoplotR::kde(x, plot = FALSE)

asVectorFunction <- function(func, onlyAllowSingleValue = TRUE) {
  if (onlyAllowSingleValue) {
    return(function(vec) list(x = vec, y = apply(as.matrix(vec), 1, func)))
  } else {
    return(function(vec) list(x = vec, y = func(vec)))
  }
}

asCumFn <- function(dat) {
  len <- length(dat)
  orderDat <- sort(dat)
  return(function(x)(sum(orderDat <= x) / len))
}

rangeToVector <- function(range, length = 400) seq(range[1], range[2], length = length)

getRangeFromDatList <- function(datList) {
  return(
    datList %>%
    lapply(function(list) range(list)) %>%
    unlist(use.names = F) %>%
    range
  )
}

getDissimilarity <- function(distanceFn) {
  return(
  function(x1, f1, x2, f2) {
    y1 <- x1 %>% f1
    y2 <- x2 %>% f2
    y1FitToX1 <- approx(x1, y1, x2, yleft = 0, yright = 1, ties = "mean")
    # print(y2 %>% length)
    # print(y1FitToX1$y %>% length)
    # plot(x1, abs(y2 - y1FitToX1$y), type="l")
    invisible(distanceFn(x2, y2, y1FitToX1$y))
  }
  )
}

getKSdistance <- function(x, y1, y2) {
  return(max(abs(y1 - y2)))
  # return(sqrt(sum((y1 - y2) ^ 2)))
}

applyListToAllCombination <- function(list, fn, isSymmetric = TRUE, isDiagonal0 = TRUE) {
  l <- length(list)
  resMatrix <- matrix(NA, nrow = l, ncol = l)
  for (i in 1:l) {
    for (j in 1:l) {
      if (!is.na(resMatrix[i, j])) next
      if (isDiagonal0 && (i == j)) {
        resMatrix[i, j] = 0;
      } else {
        resMatrix[i, j] <- fn(list[[i]], list[[j]])
        if (isSymmetric) {
          resMatrix[j, i] <- resMatrix[i, j]
        }
      }
    }
  }
  invisible(resMatrix)
}

filterList <- function(list, fn) {
  return(list[list %>% lapply(fn) %>% unlist])
}

resample <- function(vec) sample(vec, size = length(vec), replace = TRUE)

addRandomNoise <- function(vec, noise = 1) {
  return(vec + rnorm(length(vec), 0, sd = noise))
}

makeResampleDatList <- function(datList, noise = 0) {
  return(
      datList %>%
    lapply(function(vec) sample(vec, size = length(vec), replace = TRUE)) %>%
    lapply(function(vec) vec + rnorm(length(vec), 0, sd = noise))
    )
}

asContinuous <- function(vec1, vec2) {
  return(function(x) {
    if (length(vec) != length(vec)) stop('vec1 and vec2 must be same size')
    return(approx(vec1, vec2, x))
  })
}


execBootstrap <- function(dat, resampleFn, func, len = 1000) {
  res <- list()
  for (i in 1:len) {
    res[[i]] <- dat %>% resampleFn %>% func
  }
  return(res)
}


modifyMDSResult <- function(X) {
  tr <- function(m) m %>% diag %>% sum
  return(function(Yp) {
    k <- nrow(X)
    svdRes <- svd(t(X) %*% Yp)
    rT <- svdRes$v %*% t(svdRes$u)
    s <- (tr(t(X) %*% (Yp %*% rT))) / (tr(t(Yp) %*% Yp))

    # print("=----")
    # plot(X, col = 1:nrow(X), xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
    # points(Yp %*% rT,col=1:nrow(X),pch=19)
    # abline(h = 0)
    # abline(v = 0)
    # print(tr(t(rT) %*% rT))
    # print(tr(t(Yp) %*% Yp))
    # print(tr(t(Yp %*% rT) %*% (Yp %*% rT)))
    # print(tr(t(X) %*% (Yp %*% rT)))
    # print(rT)
    # print(s)
    # print(sum(colSums(((Yp %*% rT) - X)^2)))
    # print(sum(colSums((s * 0.9 * (Yp %*% rT) - X)^2)))
    # print(sum(colSums((s * (Yp %*% rT) - X)^2)))
    # print("=----")

    tt = 1 / k * colSums(X - s * Yp %*% rT)
    # stop()
    return(s * (Yp %*% rT))
    #  + (numeric(k) + 1) %*% t(tt))
  })
}

#p*q*N -> N*q*p
recombineMatrixList <- function(lis) {
  nList <- length(lis)
  nrMat <- nrow(lis[[1]])
  ncMat <- ncol(lis[[1]])
  res <- vector(mode = 'list', length = nrMat) %>% lapply(function(mat) {
    return(matrix(NA, nrow = nList, ncol = ncMat))
  })
  for (i in 1:nList) {
    for (j in 1:nrMat) {
      res[[j]][i,] <- lis[[i]][j,]
    }
  }
  return(res)
}

# plot(x, approxfun(x, ecdf(dat1Sorted)(x), method='linear')(x))

getMDS <- function(datList) {

  return(datList %>% applyListToAllCombination(function(dat1, dat2) {
    dat1Sorted <- sort(dat1)
    dat2Sorted <- sort(dat2)
    x <- datList %>% getRangeFromDatList %>% (function(p) seq(p[1], p[2], length = 400))
    getDissimilarity(getKSdistance)(
    # dat2Sorted,
    # x,
      dat2Sorted,
      ecdf(dat2Sorted),
    # approxfun(dat2Sorted, ecdf(dat2Sorted)(dat2Sorted), method='linear', yleft = 0, yright=1),
    # dat1Sorted,
    # x,
      dat1Sorted,
      ecdf(dat1Sorted)
    # approxfun(dat1Sorted, ecdf(dat1Sorted)(dat1Sorted), method='linear',yleft = 0,yright = 1)
    # ecdf(dat1Sorted)
  )
  }) %>%
  cmdscale
  # (MASS::isoMDS) %>% (function(a) a$points)
  )
}

asPairedList <- function(list1, list2, namesList1 = 'list1', namesList2 = 'list2') {
  res = vector(mode = 'list', length = length(list1))
  for (i in 1:length(list1)) {
    summary <- list()
    summary[[namesList1]] <- list1[[i]]
    summary[[namesList2]] <- list2[[i]]
    res[[i]] <- summary
  }
  return(res)
}

convertMatrixRowAsList <- function(mat) {
  res <- as.list(as.data.frame(t(mat)))
  names(res) <- NULL
  return(res)
}


getBootstrapMDS <- function(datList, ci = 0.9, Nbootstrap = 300) {
  X <- datList %>% getMDS

  bootstrapCoordinate <-
  vector(mode = 'list', length = Nbootstrap) %>% #bootstrap result container
  lapply(function(i) {
    #bootstrap process
    return(
    datList %>% #initial data set
    makeResampleDatList %>% # make bootstrap sample
    getMDS %>% # get MDS result
    (modifyMDSResult(X)) #convert the result to fit X
    )
  })

  covMat <- bootstrapCoordinate %>%
  recombineMatrixList %>%
  lapply(function(mat) cov(mat)) # make variance-covariance matrix

  ellipseCoordinate <-
  asPairedList(covMat, convertMatrixRowAsList(X), 'cov', 'mean') %>%
  lapply(function(d) {
    ellipse::ellipse(d$cov, centre = d$mean, level = ci)
  })

  return(
    list(
      X = X,
      bootstrapCoordinate = bootstrapCoordinate,
      ellipseCoordinate = ellipseCoordinate
     )
    )
}

plotBootstrapMDS <- function(resBootstrapMDS, addPoints = FALSE) {
  plot(resBootstrapMDS$X, col = 1:nrow(resBootstrapMDS$X), xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
  abline(v = 0)
  abline(h = 0)
  resBootstrapMDS$ellipseCoordinate %>% lapply(function(x) {
    lines(x)
  })
  if (addPoints) {
    resBootstrapMDS$bootstrapCoordinate %>% lapply(function(x) {
      points(x, col = 1:nrow(x), cex = .4)
    })
  }
}



# # Yp <- Yp2
# X <- datList %>% getMDS
# tr <- function(m) m %>% diag %>% sum
# k <- nrow(X)
# svdRes <- svd(t(X) %*% Yp)
# rT <- svdRes$v %*% t(svdRes$u)
# s <- (tr(t(X) %*% (Yp %*% rT))) / (tr(t(Yp) %*% Yp))
# tt = 1 / k * colSums(X - s * Yp %*% rT)
# r <- s * Yp %*% rT + (numeric(k) + 1) %*% t(tt)

# # datList %>% makeResampleDatList %>% getMDS %>%
# # %>% (modifyMDSResult(X))

# plot(X, col=1:nrow(X), xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
# points(Yp %*% rT * s, col=1:nrow(Yp), pch=19)
# abline(h = 0)
# abline(v = 0)

# # x <- datList %>% getRangeFromDatList %>% rangeToVector
# # datList %>% getMDS

# # IsoplotRMatrix <- as.matrix(IsoplotR::mds(datList)$diss)
# # cmdscale(IsoplotRMatrix)
# # IsoplotRRes <- IsoplotR::mds(datList)$points

# # MASS::isoMDS(IsoplotRMatrix)$points
# # MASS::isoMDS(dissimilarityMatrix)$points

# # plot(MASS::isoMDS(IsoplotRMatrix)$points)
# # points(MASS::isoMDS(dissimilarityMatrix)$points, pch = 10)


datList <- IsoplotR::examples$DZ %>% filterList(function(x) {
  return(TRUE)
})
# plot(1:nrow(X), 1:nrow(X), col=1:nrow(X))

datListUsed = datList %>% lapply(function(x) {
  # return (c(x,x,x))
  x
  # if(x %>% length == 100) {
  #   return (c(x,x,x,x,x))
  # } else {
  #   return(c(x,x,x,x,x,x,x,x,x,x,x,x))
  # }
  # plot(density(x))
})

datListUsed %>% lapply(function(x) {
  x %>% length %>% print
})

res <- getBootstrapMDS(datListUsed, Nbootstrap = 100)
plotBootstrapMDS(res, addPoints = TRUE)
# # %>% filterList(function(l) length(l) == 100) 
# # %>% lapply(function(x) x[c(1:floor(length(x) / 5))])

# Nbootstrap <- 400
# ci <- 0.90
