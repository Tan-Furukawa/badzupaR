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

    datList %>% applyListToAllCombination(function(dat1, dat2) {
    dat1Sorted <- sort(dat1)
    dat2Sorted <- sort(dat2)
    x <- datList %>% getRangeFromDatList %>% (function(p) seq(p[1],p[2],length=400))
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
  cmdscale -> res
  # (MASS::isoMDS) %>% (function(a) a$points)
  return(res)
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

divideDataFrameAsList <- function(df, getIndex) {
  index <- unique(getIndex(df))
  index %>% 
    lapply(function(i) {
      df[getIndex(df) == i,]
    }) -> res
  return(res)
}

reshapeList <- function (list, reduce = function(prev, curr, index) {return(curr)}) {
  res <- NULL
  for (i in 1:length(list)) {
    res <- reduce(res, list[[i]], i)
  }
  return(res)
}



#' getBootstrapMDS
#'
#' This function computes Multidimensional Scaling (MDS) and generates confidence ellipses
#' based on bootstrapped data for a list of datasets. It allows you to estimate the MDS
#' coordinates and corresponding confidence intervals.
#'
#' @param datList A list of numeric vectors representing datasets.
#' @param ci The confidence level for ellipses (default is 0.9).
#' @param Nbootstrap The number of bootstrap samples to generate (default is 100).
#' @param watcher A custom function to monitor the progress of bootstrapping.
#'
#' @return A list containing the following components:
#' \item{datList}{The input dataset list.}
#' \item{X}{The MDS coordinates of the original data.}
#' \item{bootstrapCoordinate}{A list of MDS coordinates for bootstrapped samples.}
#' \item{ellipseCoordinate}{A list of confidence ellipses for the bootstrapped data.}
#'
#' @examples
#' datList <- list(
#'   N1 = rnorm(200),
#'   N2 = rnorm(200),
#'   N3 = rnorm(200, 2, 1),
#'   N4 = rnorm(200, 2, 1),
#'   N5 = rnorm(200, 3, 1),
#'   N6 = rnorm(200, 3, 1)
#'   )
#'  res <- getBootstrapMDS(datList, Nbootstrap = 100)
#'  plotBootstrapMDS (res, addPoints=F)
#' 
#' @importFrom ellipse ellipse
#' @export
getBootstrapMDS <- function(datList, ci=0.9, Nbootstrap=100, watcher=function(progress, total) {}) {
  X <- datList %>% getMDS
  progress <- 1

  bootstrapCoordinate <-
  vector(mode = 'list', length = Nbootstrap) %>% #bootstrap result container
  lapply(function(i) {#bootstrap process
    watcher(parent.frame()$i, Nbootstrap)
    
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
    res <- ellipse::ellipse(d$cov, centre = d$mean, level = ci, npoints = 20)
    
    return(res)
  })
  
  return(
    list(
      datList = datList,
      X = X,
      bootstrapCoordinate = 
        bootstrapCoordinate %>%
        reshapeList(function(prev, curr, index) {
          return (rbind(prev, data.frame(id = numeric(nrow(curr)) + index, x = curr[,1], y = curr[,2])))
        }),
      ellipseCoordinate = ellipseCoordinate %>%
        reshapeList(function(prev, curr, index) {
          return (rbind(prev, data.frame(id = numeric(nrow(curr)) + index, x = curr[,1], y = curr[,2])))
        })
     )
    )
}

#' plotBootstrapMDS
#'
#' This function visualizes the Multidimensional Scaling (MDS) coordinates and confidence ellipses
#' generated by the `getBootstrapMDS` function. It provides a scatterplot of MDS coordinates and
#' optionally adds points from bootstrapped data along with ellipses representing confidence intervals.
#'
#' @param resBootstrapMDS A list containing MDS and ellipse coordinates obtained from `getBootstrapMDS`.
#' @param addPoints A logical value indicating whether to add bootstrapped points to the plot (default is FALSE).
#'
#' @examples
#' datList <- list(
#'   N1 = rnorm(200),
#'   N2 = rnorm(200),
#'   N3 = rnorm(200, 2, 1),
#'   N4 = rnorm(200, 2, 1),
#'   N5 = rnorm(200, 3, 1),
#'   N6 = rnorm(200, 3, 1)
#'   )
#'  res <- getBootstrapMDS(datList, Nbootstrap = 100)
#'  plotBootstrapMDS (res, addPoints=F)
plotBootstrapMDS <- function(resBootstrapMDS, addPoints = FALSE) {
  plot(resBootstrapMDS$X, col = 1:nrow(resBootstrapMDS$X), xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
  abline(v = 0)
  abline(h = 0)
  if (addPoints) {
   resBootstrapMDS$bootstrapCoordinate %>% 
      divideDataFrameAsList (function(x) x$id) %>%
      lapply(function(x) {
      points(x = x$x, y= x$y, col = 1:nrow(x), cex = .4)
    })
  }
  resBootstrapMDS$ellipseCoordinate %>% 
    divideDataFrameAsList (function(x) x$id) %>%
    lapply(function(x) {
    lines(x = x$x, y= x$y)
  })
  text(x = res$X[,1], y = res$X[,2], resBootstrapMDS$datList %>% names)
}
