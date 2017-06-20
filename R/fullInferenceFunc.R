computeWeightedPval <- function(mu, naive, contrast, cov, sampleDat,
                                smaller, larger) {
  # interpolating mean
  uweight <- (mu - smaller) / (larger - smaller)
  lweight <- 1 - uweight
  lestimate <- sampleDat$estimate[[which(sampleDat$mean == smaller)]]
  uestimate <- sampleDat$estimate[[which(sampleDat$mean == larger)]]
  estimate <- lestimate * lweight + uestimate * uweight

  # Computing weighted Sample
  if(lweight > 0.5) {
    sample <- sampleDat$samples[[which(sampleDat$mean == smaller)]]
    origDens <- mvtnorm::dmvnorm(sample, lestimate, subCov, log = TRUE)
  } else {
    sample <- sampleDat$samples[[which(sampleDat$mean == larger)]]
    origDens <- mvtnorm::dmvnorm(sample, uestimate, subCov, log = TRUE)
  }
  targetDens <- mvtnorm::dmvnorm(sample, estimate, subCov, log = TRUE)
  importanceWeights <- targetDens - origDens
  importanceWeights <- importanceWeights - max(importanceWeights)
  importanceWeights <- exp(importanceWeights)
  importanceWeights <- importanceWeights / sum(importanceWeights)

  meanSamp <- as.numeric(sample %*% contrast)
  if(naive > 0) {
    pvalue <- weighted.mean(naive < meanSamp, importanceWeights)
  } else {
    pvalue <- weighted.mean(naive > meanSamp, importanceWeights)
  }

  #print(c(mu, pvalue))
  return(pvalue)
}

computeCIquantile <- function(naive, contrast, cov, sampleDat, alpha) {
  quantiles <- c(1 - alpha / 2, alpha / 2)
  CI <- numeric(2)
  for(i in 1:2) {
    q <- quantiles[i]
    lower <- sampleDat$mean[sampleDat$pval < q][which.max(sampleDat$pval[sampleDat$pval < q])]
    upper <- sampleDat$mean[sampleDat$pval > q][which.min(sampleDat$pval[sampleDat$pval > q])]
    interval <- sort(c(lower, upper))
    c <- NULL
    try(c <- uniroot(f = function(x) computeWeightedPval(x, naive, contrast, subCov, sampleDat, lower, upper) - q,
                     interval = interval)$root)
    if(!is.null(c)) {
      CI[i] <- c
    } else {
      if(i == 1) {
        CI[i] <- max(sampleDat$mean)
      } else {
        CI[i] <- min(sampleDat$mean)
      }
    }
  }
  CI <- sort(CI)

  return(CI)
}

selectiveMRI_control <- function(barrierCoef = 0.1, stepSizeCoef = 2,
                                 stepRate = 0.6, trimSample = 100,
                                 delay = 10, optimSteps = 500,
                                 CIstepSize = 1, maxCItries = NULL,
                                 nSamples = NULL) {
  control <- list(barrierCoef = barrierCoef,
                  stepSizeCoef = stepSizeCoef,
                  stepRate = stepRate,
                  trimSample = trimSample,
                  delay = delay,
                  optimSteps = optimSteps,
                  CIstepSize = CIstepSize,
                  maxCItries = maxCItries,
                  nSamples = nSamples)
  class(control) <- "fmriControl"
  return(control)
}

selectiveMRI <- function(y, cov, threshold,
                         contrast = NULL,
                         coordinates = NULL,
                         selected = NULL,
                         tykohonovSlack = 1,
                         CIalpha = 0.05,
                         probMethod = c("all", "selected", "onesided"),
                         imputeBoundary = c("none", "mean", "neighbors"),
                         control = selectiveMRI_control()) {
  # Preparations ------------------
  CIstepSize <- control$CIstepSize
  maxCItries <- control$maxCItries
  if(is.null(maxCItries)) {
    maxCItries <- 10 / CIstepSize
  }
  barrierCoef <- control$barrierCoef
  stepSizeCoef <- control$stepSizeCoef
  stepRate <- control$stepRate
  trimSample <- control$trimSample
  delay <- control$delay
  optimSteps <- control$optimSteps
  nSamples <- control$nSamples
  if(is.null(nSamples)) {
    nSamples <- 1 / CIalpha * 20
  }
  if(is.null(contrast)) {
    contrast <- rep(0, nrow(coordinates))
    contrast[selected] <- rep(1 / sum(selected), sum(selected))
  }

  # Computing MLE -------------------------
  print("Computing MLE!")
  mle <- optimizeSelected(y, cov, threshold,
                          selected = selected,
                          projected = NULL,
                          stepRate = stepRate,
                          coordinates = coordinates,
                          tykohonovParam = NULL,
                          tykohonovSlack = tykohonovSlack,
                          stepSizeCoef = stepSizeCoef,
                          delay = delay,
                          assumeConvergence = optimSteps,
                          trimSample = trimSample,
                          maxiter = optimSteps + 1,
                          probMethod = probMethod,
                          imputeBoundary = imputeBoundary)

  conditional <- sum(mle$conditional * contrast)
  naive <- sum(y * contrast)

  # p-value at MLE
  samp <- mle$sample
  mlemeans <- as.numeric(samp %*% contrast)
  if(naive > 0) {
    mlepval <- mean(mlemeans > naive)
  } else {
    mlepval <- mean(mlemeans < naive)
  }

  # Computing p-value ----------------------------------------
  print("Computing p-value!")
  nullfit <- optimizeSelected(y, cov, threshold,
                              projected = 0,
                              selected = selected,
                              stepRate = stepRate,
                              coordinates = coordinates,
                              tykohonovParam = NULL,
                              tykohonovSlack = 0.0001,
                              stepSizeCoef = 0,
                              delay = 1,
                              assumeConvergence = 1,
                              trimSample = trimSample,
                              maxiter = nSamples*2 + 1,
                              probMethod = probMethod,
                              init = rep(0, length(y)),
                              imputeBoundary = "neighbors")
  samp <- nullfit$sample
  nullmeans <- as.numeric(samp %*% contrast)
  if(naive > 0) {
    pvalue <- mean(naive < nullmeans)
  } else {
    pvalue <- mean(naive > nullmeans)
  }

  # Setting up `sample dataset` ---------------------
  profSamples <- data.frame(mean = c(0, conditional),
                            estimate = I(list(rep(0, length(y)), mle$conditional)),
                            samples = I(list(nullfit$sample, mle$sample)),
                            pval = c(pvalue, mlepval))

  print("Finding lower CI bound!")
  # Find lower boundary for search -----------------
  if(pvalue > CIalpha) {
    projMean <- 0
    direction <- -1
    test <- pvalue
    init <- rep(0, length(y))
  } else {
    projMean <- conditional
    direction <- -1
    test <- mlepval
    init <- mle$conditional
  }
  test <- pvalue
  meansd <- as.numeric(t(contrast) %*%  subCov %*% contrast)
  try <- 1
  while(test > CIalpha / 2 & try <= maxCItries) {
    projMean <- projMean + direction * meansd * sign(naive)
    try(profile <- optimizeSelected(y, cov, threshold,
                                    selected = selected,
                                    projected = projMean,
                                    stepRate = stepRate,
                                    coordinates = coordinates,
                                    tykohonovSlack = tykohonovSlack,
                                    stepSizeCoef = stepSizeCoef,
                                    delay = delay,
                                    assumeConvergence = optimSteps / 2,
                                    trimSample = trimSample,
                                    maxiter = optimSteps / 2 + nSamples,
                                    probMethod = probMethod,
                                    init = NULL,
                                    imputeBoundary = imputeBoundary))
    samp <- profile$sample
    profmeans <- as.numeric(samp %*% contrast)
    if(naive > 0) {
      test <- mean(profmeans > naive)
    } else {
      test <- mean(profmeans < naive)
    }
    profSamples <- rbind(profSamples, data.frame(mean = projMean,
                                                 estimate = I(list(profile$conditional)),
                                                 samples = I(list(profile$sample)),
                                                 pval = test))
    try <- try + 1
    print(c(-1, try, projMean, test))
  }

  print("Finding upper CI bound!")
  # Finding upper bound for CI search --------------------
  try(profile <- optimizeSelected(y, cov, threshold,
                                  selected = selected,
                                  projected = naive,
                                  stepRate = stepRate,
                                  coordinates = coordinates,
                                  tykohonovSlack = tykohonovSlack,
                                  stepSizeCoef = stepSizeCoef,
                                  delay = delay,
                                  assumeConvergence = optimSteps / 2,
                                  trimSample = trimSample,
                                  maxiter = optimSteps / 2 + nSamples,
                                  probMethod = probMethod,
                                  init = NULL,
                                  imputeBoundary = imputeBoundary))
  samp <- profile$sample
  profmeans <- as.numeric(samp %*% contrast)
  if(naive > 0) {
    test <- mean(profmeans > naive)
  } else {
    test <- mean(profmeans < naive)
  }

  profSamples <- rbind(profSamples, data.frame(mean = naive,
                                               estimate = I(list(profile$conditional)),
                                               samples = I(list(profile$sample)),
                                               pval = test))
  projMean <- naive
  if(test > 1 - alpha / 2) {
    direction <- -1
    conditionf <- function(x) x > 1 - CIalpha / 2
  } else {
    direction <- 1
    conditionf <- function(x) x < 1 - CIalpha / 2
  }

  try <- 1
  while(conditionf(test) & try <= maxCItries) {
    projMean <- projMean + direction * meansd * sign(naive)
    try(profile <- optimizeSelected(y, cov, threshold,
                                    selected = selected,
                                    projected = projMean,
                                    stepRate = stepRate,
                                    coordinates = coordinates,
                                    tykohonovSlack = tykohonovSlack,
                                    stepSizeCoef = stepSizeCoef,
                                    delay = delay,
                                    assumeConvergence = optimSteps / 2,
                                    trimSample = trimSample,
                                    maxiter = optimSteps / 2 + nSamples,
                                    probMethod = probMethod,
                                    init = NULL,
                                    imputeBoundary = imputeBoundary))
    samp <- profile$sample
    profmeans <- as.numeric(samp %*% contrast)
    if(naive > 0) {
      test <- mean(profmeans > naive)
    } else {
      test <- mean(profmeans < naive)
    }
    profSamples <- rbind(profSamples, data.frame(mean = projMean,
                                                 estimate = I(list(profile$conditional)),
                                                 samples = I(list(profile$sample)),
                                                 pval = test))
    try <- try + 1
    print(c(-1, try, projMean, test))
  }

  print("Done!")
  CI <- computeCIquantile(naive, contrast, subCov, profSamples, CIalpha)
  return(list(conditional = conditional,
              CI = CI,
              sampleDat = profSamples))
}
