run.sim <- function(config, noise_type ="sim", noise_dat = NULL) {

  snr <- config[["snr"]]
  spread <- config[["spread"]]
  BHlevel <- config[["BHlevel"]]
  replications <- config[["replications"]]

  if(noise_type == "sim" ) {
    rho <- config[["rho"]]
    grp_size <- -1
  } else {
    grp_size <- config[["grp_size"]]
    rho <- -1
  }

  I <- 11
  J <- 11
  K <- 9
  dims <- c(I, J, K)

  coordinates <- expand.grid(i = 1:I, j = 1:J, k = 1:K)

  if (noise_type == "sim") {
    covariance <- rho^as.matrix(dist(coordinates[, 1:3], method = "euclidean",
                                   diag = TRUE, upper = TRUE))
    covEigen <- eigen(covariance)
    sqrtCov <- covEigen$vectors %*% diag(sqrt(covEigen$values)) %*% t(covEigen$vectors)

  }
  else if (noise_type == "fmri"){
    stopifnot(all.equal(dim(noise_dat)[1:3],dims))
    pop_size <- dim(noise_dat)[4]
    dim(noise_dat) = c(prod(dims), pop_size)
    alpha <- 1/pop_size
    sds <- apply(noise_dat, 1, sd)
    # covariance of a single sample
    sample_covariance <- (1-alpha)*cov(t(noise_dat))+alpha*diag(sds^2)
    # covariance of group averages
    covariance <- sample_covariance*(1/grp_size*2)
    dim(noise_dat) <- c(dims[1:3],pop_size)
    rm(sample_covariance)
  }


  simresults <- list()
  simcover <- matrix(nrow = replications, ncol = 2)
  colnames(simcover) <- c("cover", "power")
  for(rep in 1:replications) {
    selected <- FALSE
    maxsize <- 0
    while(is.na(maxsize) | maxsize < 2) {
      if (noise_type == "sim") {
        coordinates <- generateArrayData3D(dims, sqrtCov, snr, spread)
      }
      else if (noise_type == "fmri"){
        coordinates <- residualData3D(noise_dat, grp_size, snr, spread)
      }
      coordinates$zval <- coordinates$observed / sqrt(diag(covariance))
      coordinates$pval <- 2 * pnorm(-abs(coordinates$zval))
      coordinates$qval <- p.adjust(coordinates$pval, method = "bonferroni")
      coordinates$qval <- coordinates$pval  # FOR NO CORRECTION!!!!
      coordinates$selected <- coordinates$qval < BHlevel
      selected <- coordinates$selected
      if(sum(selected) >= 2) {
        clusters <- findClusters(coordinates)
        maxsize <- max(sapply(clusters, function(x) sum(x$selected)))
      }
    }

    print(c(rep = rep, rho = rho, grp_size = grp_size, BHlevel = BHlevel, snr = snr, spread = spread))
    print(c(nselected = sum(selected)))
    sizes <- sapply(clusters, nrow)
    # threshold <- qnorm(BHlevel * sum(coordinates$selected) / nrow(coordinates) / 2,
    #                    lower.tail = FALSE) ## This is BH threshold
    #threshold <- qnorm(1 - BHlevel / (I * J * K) / 2)
    threshold <- abs(qnorm(BHlevel / 2))

    results <- list()
    iterPower <- 0
    iterCover <- 0
    weights <- 0
    slot <- 1
    mse <- 0
    msecount <- 0
    for(m in 1:length(clusters)) {
      cluster <- clusters[[m]]
      cluster <- subset(cluster, !is.na(cluster$selected))
      selected <- coordinates$selected[cluster$row]
      observed <- coordinates$observed[cluster$row]
      signal <- coordinates$signal[cluster$row]
      subCov <- covariance[cluster$row, cluster$row, drop = FALSE]

      print(c(round(m / length(clusters), 2), nrow(cluster)))

      try(mle <- optimizeSelected(observed, subCov, threshold,
                                  selected = selected,
                                  projected = NULL,
                                  stepRate = 0.6,
                                  coordinates = cluster[, 1:3],
                                  tykohonovParam = NULL,
                                  tykohonovSlack = 2,
                                  stepSizeCoef = 2,
                                  delay = 10,
                                  assumeConvergence = 1000,
                                  trimSample = 200,
                                  maxiter = 1200,
                                  probMethod = "selected",
                                  init = observed,
                                  imputeBoundary = "neighbors"))

      conditional <- mean(mle$conditional[selected])
      naive <- mean(observed[selected])
      true <- mean(signal[selected])
      mse <- mse * msecount / (msecount + 1) + c(naive = (naive - true)^2, conditional = (conditional - true)^2) / (msecount + 1)
      msecount <- msecount + 1
      for(methodind in 1:2) {
        if(methodind == 1) {
          method <- "onesided"
        } else if(methodind == 2) {
          method <- "selected"
        }
        results[[slot]] <- list()

        # Computing the p-value based on samples from the null
        nullfit <- NULL
        try(nullfit <- optimizeSelected(observed, subCov, threshold,
                                       projected = 0,
                                       selected = selected,
                                       stepRate = 0.6,
                                       coordinates = cluster[, 1:3],
                                       tykohonovParam = NULL,
                                       tykohonovSlack = 0.0001,
                                       stepSizeCoef = 0,
                                       delay = 1,
                                       assumeConvergence = 1,
                                       trimSample = 200,
                                       maxiter = 1500,
                                       probMethod = method,
                                       init = rep(0, length(observed)),
                                       imputeBoundary = "neighbors"))
        if(is.null(nullfit)) next
        naive <- mean(observed[selected])
        samp <- nullfit$sample
        nullmeans <- NULL
        try(nullmeans <- rowMeans(samp[, selected, drop = FALSE]))
        if(is.null(nullmeans)) next
        # p-value based on one sided test
        if(naive > 0) {
          pvalue <- mean(naive < nullmeans)
        } else {
          pvalue <- mean(naive > nullmeans)
        }

        # Projecting to the truth, rejecting the test here implies
        # that the CI doesn't cover the truth.
        try(profile <- optimizeSelected(observed, subCov, threshold,
                                        selected = selected,
                                        projected = mean(signal[selected]),
                                        stepRate = 0.6,
                                        coordinates = cluster[, 1:3],
                                        tykohonovParam = NULL,
                                        tykohonovSlack = 2,
                                        stepSizeCoef = 2,
                                        delay = 10,
                                        assumeConvergence = 750,
                                        trimSample = 200,
                                        maxiter = 2000,
                                        probMethod = method,
                                        init = observed,
                                        imputeBoundary = "neighbors"))

        samp <- profile$sample
        profMeans <- NULL
        try(profMeans <- rowMeans(samp[, selected, drop = FALSE]))
        if(is.null(profMeans)) next
        # Two sided CI
        profPval <- 2 * min(mean(naive < profMeans), mean(naive > profMeans))

        true <- mean(signal[selected])
        profResult <- c(true = mean(signal[selected]), profPval = profPval, pvalue = pvalue)
        results[[slot]][[1]] <- c(snr = snr, spread = spread, method = methodind,
                               rho = rho, BHlevel = BHlevel,
                               size = sum(selected), true = true)
        results[[slot]][[2]] <- profResult
        results[[slot]][[3]] <- c(conditional = conditional, naive = naive, true = true)

        print(profResult)
        print(mse)

        weight <- 1
        iterPower <- iterPower + weight * (pvalue < 0.05)
        iterCover <- iterCover + weight * (profPval > 0.05)
        weights <- weights + weight

        # plot(density((rowMeans(profile$sample[, selected, drop = FALSE]))), xlim = c(-4.5, 4.5))
        # abline(v = mean(profile$conditional[selected]))
        # abline(v = mean(result$conditional[selected]), col = "blue")
        # abline(v = naive, col = "red")
        # lines(density((rowMeans(result$sample[, selected]))), col = "blue")
        # abline(v = c(threshold, - threshold), col = "grey")
        slot <- slot + 1
      }
    }

    simcover[rep, 1] <- sum(iterCover) / sum(weights)
    simcover[rep, 2] <- sum(iterPower) / sum(weights)
    print(c(cover = iterCover, power = iterPower) / weights)
    print(colMeans(simcover[1:rep, , drop = FALSE]))
    simresults[[rep]] <- results
  }

  return(simresults)
}

configurations <- expand.grid(snr = c(4, 3, 2, 1, 0),
                              spread = c(2),
                              rho = c(0.5, 0.75),
                              BHlevel = c(0.001, 0.01),
                              replications = 2,
                              grp_size = c(4,16))

system.time(simResults <- apply(configurations, 1, run.sim, noise_type ="sim"))
save(simResults, file = "simulations/results/May 18 sim.Robj")

load('/Users/yuvalb/Dropbox/SelectiveFmri/brain_data_4mm_Cambridge.rda')
dat = brain_data$t_cube[1:11,1:11,1:9,]

system.time(simResults <- apply(configurations, 1, run.sim, noise_type ="fmri",dat))
save(simResults, file = "simulations/results/May 18 fmri.Robj")




