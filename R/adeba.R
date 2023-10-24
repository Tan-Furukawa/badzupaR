# Rcpp::sourceCpp('src/make_pilot.cpp')
# Rcpp::sourceCpp('src/CppLogLikelihood.cpp')
# library("ramify")

get_bandwidth = function(alpha, beta, pilot) {
    return(alpha / pilot ^ beta)
}

roop_likelihoood_by_alpha = function(beta, vec_alpha, pilot, dat) {
    dat <- sort(dat)
    discre_size <- length(vec_alpha)
    n <- length(dat)
    likelihood_by_alpha <- numeric(discre_size)
    for (i in 1:length(vec_alpha)) {
        bandwidth <- get_bandwidth(vec_alpha[i], beta, pilot)
        likelihood_by_alpha[i] <-
            CppLoglikelihood(dat, bandwidth, as.integer(n))
    }
    return(likelihood_by_alpha)
}

argmax_logalpha_in_likelihood = function(range_log_alpha,
                                            pilot,
                                            dat,
                                            discre_size,
                                            beta) {
    vec_logalpha <- seq(range_log_alpha[1], range_log_alpha[2], length = discre_size + 1)[1:discre_size]
    vec_alpha <- 10 ^ vec_logalpha
    likelihood_by_alpha <-
        roop_likelihoood_by_alpha(beta, vec_alpha, pilot, dat)
    
    argmax <- ramify::argmax(t(as.matrix(likelihood_by_alpha)))
    begin <- NA
    end <- NA
    l_alpha <- length(likelihood_by_alpha)
    if (argmax == 1) {
        begin <- 1
        end <- 1
        warning("range of alpha is not enough")
    } else if (argmax == l_alpha) {
        begin <- l_alpha
        end <- l_alpha
        warning("range of alpha is not enough")
    } else {
        begin <- argmax - 1
        end <- argmax + 1
    }
    
    return(
        list(
            argmax = argmax,
            max_alpha = vec_logalpha[argmax],
            range_log_alpha = c(vec_logalpha[begin], vec_logalpha[end]),
            likelihood = likelihood_by_alpha[argmax]
        )
    )
}

get_adeba_bandwidth = function(dat,
                                range_log_alpha,
                                vec_beta,
                                discre_size1,
                                discre_size2,
                                alpha0,
                                beta0) {
    pilot0 <- numeric(length(dat)) + alpha0
    res1 <-
        argmax_logalpha_in_likelihood(range_log_alpha, pilot0, dat, discre_size1, beta0)
    
    res2 <-
        argmax_logalpha_in_likelihood(res1$range_log_alpha, pilot0, dat, discre_size2, beta0)
    argmax_alpha_result <- res2$argmax
    
    n <- length(dat)
    likelihood_by_beta <-
        numeric(length(vec_beta)) - Inf
    pilot1_prev <- numeric(n)
    alpha1 <- 10 ^ res2$max_alpha
    
    for (i in 1:length(vec_beta)) {
        bandwidth <- get_bandwidth(alpha1, vec_beta[i], pilot0)
        pilot1 <-
            make_pilot(dat, bandwidth, as.integer(n))
        
        res1 <-
            argmax_logalpha_in_likelihood(range_log_alpha, pilot1, dat, discre_size1, vec_beta[i])
        res2 <-
            argmax_logalpha_in_likelihood(res1$range_log_alpha, pilot1, dat, discre_size2, vec_beta[i])
        
        likelihood <- res2$likelihood
        likelihood_by_beta[i] <- likelihood
        alpha1 <- 10 ^ res2$max_alpha
        
        argmax_beta <-
            ramify::argmax(t(as.matrix(likelihood_by_beta)))
        max_likelihood <-
            likelihood_by_beta[argmax_beta]
        beta1 <-
            vec_beta[ramify::argmax(t(as.matrix(likelihood_by_beta)))]
        
        if (max_likelihood != likelihood) {
            pilot1 <- pilot1_prev
            alpha1 <- alpha1_prev
            break
        }
        
        pilot1_prev <- pilot1
        alpha1_prev <- alpha1
    }
    
    bandwidth_result <-
        get_bandwidth(alpha1, beta1, pilot1)
    return(bandwidth_result)
}

kernel_density = function(dat,
                            bandwith,
                            kernelFn,
                            x = NULL,
                            out = .1,
                            m = 400) {
    n <- length(dat)
    delta <- max(dat) - min(dat)
    # if (is.null(x)) {
    #     x <- seq(min(dat) - delta * out, max(dat) + delta * out, length = m)
    # }
    y <- numeric(length(x))
    for (i in 1:length(dat)) {
        if (bandwith[i] == Inf) {
            next
        }
        y <-
            y + 1 / n / bandwith[i] * kernelFn((dat[i] - x) / bandwith[i])
    }
    return(list(x = x, y = y))
}

triweight = function(u) {
    res <- numeric(length(u))
    LogiAbsU <- abs(u) <= 1
    res[LogiAbsU] <-
        35 / 32 * (1 - u[LogiAbsU] ^ 2) ^ 3
    return(res)
}

#' Calculate Probability Density Using ADEBA Algorithm
#'
#' This function calculates the probability density of input data using the ADEBA
#' (Automatic DEnsity Bandwidth Adjustment) algorithm. It employs a combination of
#' alpha and beta parameters to adjust the bandwidth for kernel density estimation.
#'
#' @param x A numeric vector representing the input data.
#' @param range_log_alpha A numeric vector of length 2 indicating the search range
#'   for the parameter alpha. Increase the range if computation fails.
#' @param vec_beta A numeric vector containing candidate values for the parameter beta.
#' @param discre_size1 The size of discretization for alpha values.
#' @param discre_size2 Two values from discre_size1 are further discretized using discre_size2
#'   around the maximum likelihood estimate of alpha to improve precision.
#' @param alpha0 The initial value for the alpha parameter.
#' @param beta0 The initial value for the beta parameter.
#' @param from The lower bound of the density estimation range.
#' @param to The upper bound of the density estimation range.
#' @param out The extra padding around the data range for drawing, if from and to are missing.
#' @param n The number of points to evaluate the density.
#' @return A "density" object containing x (density values on the x-axis) and y (normalized density values).
#' @examples
#' data <- rnorm(100)
#' density_result <- density.adeba(data)
#' plot(density_result, main = "Density Estimation")
#'
#' @importFrom Rcpp evalCpp
#' @export
density.adeba = function(x,
                            range_log_alpha = c(-6.0, 1.0),
                            vec_beta = c(0.0, 0.5, 1.0, 1.5, 2.0),
                            discre_size1 = 14,
                            discre_size2 = 7,
                            alpha0 = 1.0,
                            beta0 = 0.01,
                            from = NA,
                            to = NA,
                            out = .2,
                            n = 400) {
    dat <- x
    dat <- sort(dat)
    mind <- min(dat)
    maxd <- max(dat)
    
    if (is.na(to)) {
        to <- max(dat) + (max(dat) - min(dat)) * out
    } else {
        to <- to
    }
    
    if (is.na(from)) {
        from <- min(dat) - (max(dat) - min(dat)) * out
    } else {
        from <- from
    }
    
    n_from <- (from - mind) / (maxd - mind)
    n_to <- (to - mind) / (maxd - mind)
    n_dat <- (dat - mind) / (maxd - mind)

    n_x <- seq(n_from, n_to, length = n)
    
    bandwith <- get_adeba_bandwidth(
        n_dat,
        range_log_alpha,
        vec_beta,
        discre_size1,
        discre_size2,
        alpha0,
        beta0
    )
    
    normalize_dens <- kernel_density(n_dat,
                            bandwith,
                            triweight,
                            x = n_x,
                            out = out,
                            m = n)
                            
    nx <- normalize_dens$x
    ny <- normalize_dens$y
    
    x <- nx * (maxd - mind) + mind
    delta <- (x[length(x)] - x[1])/length(x)
    sumP <- sum(ny)
    y <- ny / (sumP * delta)

    name <- deparse(dat)
    structure(list(x = x, y = y, n = length(x), 
    data = dat, call = match.call(),
    data.name = name, has.na = FALSE), 
            class = "density")
}


#' Calculate Probability Density Using ADEBA Algorithm
#' same as density.adeba
#' @export
adeba_denstiy <- function(...) {
    density.adeba (...)
}
