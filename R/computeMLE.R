# y = vector of observations
# cov = covariance estimate
# threshold = a single threshold for all observations, selection event is y > threshold | y < threshold.
# projected = Should the mean be constrained to a value, if projected is a scalar then
#             the mle is estimated constrained to \bar{mu} == projected.
# quadraticSlack = For projected gradient, the mean vector is constrained to a ball.
#                  This determines the size of the ball.
# barrierCoef = For non-constrained MLE the MLE is estimated with a barrier function preventing
#               the coordinates of mu from changing signs, this coefficient is related to the barrier function.
# stepSizeCoef = Scaling for stochastic gradient.
# stepRate = Stochastic gradient is scaled by iteration^{-stepRate}
# trimSample = How many Gibbs cycles to go through between samples.
# lambdaStart = For projected gradient method, what should the initial value for lambda
#               be. Should be set to a high value if mean is constrained to a small ball.
# maxiter = Number of stochastic gradient steps. 1000 is usually more than enough. For some projected
#           problems a large number of steps is required to properly tune lambda.

# Outputs:
# sample - a matrix of samples
# estimates - the optimization path
# conditional - the conditional estimate
# coordinateCI - confidence intervals for the coordinates
# meanCI - CI for the mean
optimizeSelected <- function(y, cov, threshold,
                             selected = NULL,
                             projected = NULL,
                             quadraticSlack = 0.1,
                             barrierCoef = 0.5,
                             stepSizeCoef = 0.25,
                             stepRate = 0.65,
                             trimSample = 40,
                             lambdaStart = 0,
                             delay = 100,
                             maxiter = 4000,
                             assumeConvergence = 2000,
                             CIalpha = 0.05,
                             init = NULL,
                             probMethod = c("all", "selected", "oneside"),
                             imputeBoundary = TRUE,
                             neighbors = NULL) {
  # Basic checks and preliminaries ---------
  if(length(probMethod) > 1) probMethod <- probMethod[1]

  if (is.null(selected)){
    selected <- abs(y) > threshold
  }
  nselected <- sum(selected)
  p <- length(y)
  s <- sum(selected)
  if(!all(sign(y[selected]) == sign(y[selected][1]))) stop("Signs of active set must be identical")

  vars <- diag(cov)
  cov <- cov2cor(cov)
  y <- y / sqrt(vars)
  if(!is.null(init)) {
    if(length(init) != length(y)) stop("If init is not NULL, then its length must match the length of y.")
    mu <- init
  } else {
    mu <- y
    mu[!selected] <- 0
  }
  threshold <- threshold / sqrt(vars)

  sds <- sqrt(diag(cov))
  invcov <- solve(cov)
  suffStat <- as.numeric(invcov %*% y)
  a <- -threshold
  b <- threshold
  signs <- sign(y)

  # Setting up projection -----------------
  # If a projected gradient method is used then initalization must be from
  # a parameter which satisfies the constraint. Also, in such a case we
  # use a quadratic barrier function. ball is the size radius in which mu
  # is allowed to be
  if(!is.null(projected)) {
    mu <- mu * sqrt(vars)
    mu[selected] <- rep(projected, sum(selected))
    mu <- mu / sqrt(vars)
    ball <- sum(selected * mean(mu[selected])^2) + quadraticSlack
  }

  estimates <- matrix(nrow = maxiter, ncol = p)
  sampleMat <- matrix(nrow = maxiter - 1, ncol = length(y))
  estimates[1, ] <- mu
  # initial values for positive/negative Gibbs samplers
  posSamp <- abs(y)
  negSamp <- -abs(y)
  if(is.null(lambdaStart)) {
    lambda <- 0
  } else {
    lambda <- lambdaStart
  }

  restarts <- 0
  for(i in 2:maxiter) {
    # Every now and then recompute probability to be positive/negative
    if((i == 2 | i < (100 + delay) & (i %% 4 == 0)) | (i %% 50 == 0 & i <= assumeConvergence)) {
      a[selected] <- threshold[selected]
      b[selected] <- Inf
      a[!selected] <- -Inf
      b[!selected] <- threshold[!selected]
      if(probMethod == "selected") {
        posProb <- mvtnorm::pmvnorm(a[selected], b[selected],
                                    mean = mu[selected],
                                    sigma = cov[selected, selected])[[1]]
      } else if(probMethod == "all") {
        posProb <- mvtnorm::pmvnorm(a, b, mu, sigma = cov)[[1]]
      } else if(probMethod == "onesided") {
        posProb <- as.numeric(sign(y[selected][1]) == 1)
      }
      a[selected] <- -Inf
      b[selected] <- -threshold[selected]
      a[!selected] <- -threshold[!selected]
      b[!selected] <- Inf
      if(probMethod == "selected") {
        negProb <- mvtnorm::pmvnorm(a[selected], b[selected],
                                    mean = mu[selected],
                                    sigma = cov[selected, selected])[[1]]
      } else  if(probMethod == "all") {
        negProb <- mvtnorm::pmvnorm(a, b, mean = mu, sigma = cov)[[1]]
      } else if(probMethod == "onesided") {
        negProb <- as.numeric(sign(y[selected][1]) == -1)
      }
      posProb <- posProb / (posProb + negProb)
    }

    # Sampling sign followed by sampling truncated normal rv
    # we used pass by reference to modify posSamp/negSamp in place.
    sampSign <- - 1 + 2 * rbinom(1, 1, posProb)
    if(sampSign == 1) {
      a[selected] <- threshold[selected]
      b[selected] <- Inf
      a[!selected] <- -Inf
      b[!selected] <- threshold[!selected]
      newsamp <- sampleTruncNorm(posSamp, a, b, mu, invcov, trimSample)
      attempts <- 0
      while(any(is.nan(newsamp))) {
        restarts <- restarts + 1
        attempts <- attempts + 1
        newsamp <- sampleTruncNorm(y, a, b, mu, invcov, 200)
        if(attempts > 100) stop("Can't sample truncated normal samples!")
      }
      posSamp <- newsamp
      samp <- newsamp
    } else {
      a[selected] <- -Inf
      b[selected] <- -threshold[selected]
      a[!selected] <- -threshold[!selected]
      b[!selected] <- Inf
      newsamp <- sampleTruncNorm(negSamp, a, b, mu, invcov, trimSample)
      attempts <- 0
      while(any(is.nan(newsamp))) {
        attempts <- attempts + 1
        restarts <- restarts + 1
        newsamp <- sampleTruncNorm(y, a, b, mu, invcov, 200)
        if(attempts > 100) stop("Can't sample truncated normal samples!")
      }
      negSamp <- newsamp
      samp <- newsamp
    }

    sampleMat[i - 1, ] <- samp
    condExp <- as.numeric(invcov %*% samp)

    # Computing barrier part of gradient
    if(is.null(projected)) {
      barrierGrad <- (sign(mu) / abs(mean(mu[selected]))) * (barrierCoef / min(max(i - delay, 1), assumeConvergence - delay)) / sum(selected)
    } else {
      lambda <- max(lambda + (sum(mu[selected]^2) - ball) * 4, 0)
      barrierGrad <- - mu * lambda * stepSizeCoef
    }

    # Computing gradient and projecting if necessary
    # The projection of the gradient is simply setting its mean to zero
    gradient <- (suffStat - condExp + barrierGrad) / max(i - delay, 1) * stepSizeCoef
    gradient[!selected] <- 0
    gradsign <- sign(gradient)
    gradient <- pmin(abs(gradient), 0.05) * gradsign
    if(!is.null(projected)) {
      gradient <- gradient * sqrt(vars)
      gradMean <- mean(gradient[selected])
      gradient[selected] <- gradient[selected] - gradMean
      gradient <- gradient / sqrt(vars)
    }

    # Updating estimate. The error thing is to make sure we didn't accidently
    # cross the barrier. There might be a better way to do this.
    mu <- mu + gradient

    #### EXPERIMENTAL #######
    if(imputeBoundary) {
      if(is.null(neighbors)) {
        mu[!selected] <- mean(mu[selected])
      } else {
        mu[neighbors[, 1]] <- mu[neighbors[, 2]]
      }
    }
    #########################

    if(is.null(projected)) {
      mu <- pmin(abs(mu), abs(y)) * sign(y)
    }
    estimates[i,] <- mu

    # progress...
    #if((i %% 100) == 0) cat(i, " ", round(mean(mu[selected]), 3), " ")
  }
  #cat("\n")

  if(restarts > 0) {
    warning(paste("Chain restarted", restarts, "times!"))
  }

  # Computations for confidence intervals ---------------------------
  forQuants <- sampleMat[assumeConvergence:(maxiter - 1), ]
  forQuants <- forQuants %*% invcov
  condVar <- var(forQuants)
  centers <- colMeans(forQuants)
  barrierDeriv <- 0#(signs / abs(mu)) * (barrierCoef / (assumeConvergence - delay))
  forQuants <- t(t(forQuants) - centers + barrierDeriv)
  sandwich <- diag(ncol(forQuants))
  sandwich[selected, ] <- condVar[selected, ]
  B <- t(forQuants) %*% forQuants / nrow(forQuants)
  S <- solve(sandwich)
  S <- t(S) %*% B %*% S
  forQuants <- forQuants %*% t(solve(sandwich))
  ciQuantiles <- apply(forQuants, 2, function(x) quantile(x, c(1 - CIalpha / 2, CIalpha / 2)))
  ciQuantiles <- t(t(ciQuantiles) * sqrt(vars))

  forQuants <- rowMeans(t(t(forQuants) * sqrt(vars))[, selected, drop = FALSE])
  meanquantiles <- quantile(forQuants, c(1 - CIalpha / 2, CIalpha / 2))

  # Unnormalizing estimates and samples --------------------------
  for(i in 1:ncol(sampleMat)) {
    sampleMat[, i] <- sampleMat[, i] * sqrt(vars[i])
    estimates[, i] <- estimates[, i] * sqrt(vars[i])
  }

  # Computing estimate and CIs -----------------------------
  conditional <- colMeans(estimates[floor(maxiter * 0.8):maxiter, ])
  meanCI <- mean(conditional[selected]) - meanquantiles

  CI <- matrix(nrow = length(conditional), ncol = 2)
  for(i in 1:ncol(ciQuantiles)) {
    CI[i, ] <- conditional[i] - ciQuantiles[, i]
  }
  CI <- CI[selected, , drop = FALSE]

  return(list(sample = sampleMat,
              estimates = estimates,
              conditional = conditional,
              coordinateCI = CI,
              meanCI = meanCI))
}
