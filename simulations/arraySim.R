run.sim <- function(config) {
  snr <- config[[1]]
  rho <- config[[2]]
  BHlevel <- config[[3]]
  replications <- config[[4]]

  I <- 11
  J <- 11
  K <- 9
  dims <- c(I, J, K)

  coordinates <- expand.grid(i = 1:I, j = 1:J, k = 1:K)
  covariance <- rho^as.matrix(dist(coordinates[, 1:3], method = "euclidean",
                                   diag = TRUE, upper = TRUE))
  covEigen <- eigen(covariance)
  sqrtCov <- covEigen$vectors %*% diag(sqrt(covEigen$values)) %*% t(covEigen$vectors)

  simresults <- list()
  simcover <- matrix(nrow = replications, ncol = 3)
  colnames(simcover) <- c("naive", "mle", "profile")
  for(rep in 1:replications) {
    selected <- FALSE
    while(sum(selected) < 5) {
      coordinates <- generateArrayData3D(dims, sqrtCov, snr)
      coordinates$zval <- coordinates$observed / sqrt(diag(covariance))
      coordinates$pval <- 2 * pnorm(-abs(coordinates$zval))
      coordinates$qval <- p.adjust(coordinates$pval, method = "BH")
      coordinates$selected <- coordinates$qval < BHlevel
      selected <- coordinates$selected
    }

    print(c(rep = rep, rho = rho, BHlevel = BHlevel, snr = snr))
    print(c(nselected = sum(selected)))
    clusters <- findClusters(coordinates)
    sizes <- sapply(clusters, nrow)
    threshold <- qnorm(BHlevel * sum(coordinates$selected) / nrow(coordinates) / 2,
                       lower.tail = FALSE)

    results <- list()
    iterCover <- 0
    naiveCover <- 0
    profCover <- 0
    weights <- 0
    for(m in 1:length(clusters)) {
      results[[m]] <- list()
      cluster <- clusters[[m]]
      cluster <- subset(cluster, !is.na(cluster$selected))
      if(nrow(cluster) == 1) next
      print(c(round(m / length(clusters), 2), nrow(cluster)))
      subCov <- covariance[cluster$row, cluster$row, drop = FALSE]

      observed <- coordinates$observed[cluster$row]
      result <- NULL
      try(result <- optimizeSelected(observed, subCov, threshold,
                                 stepRate = 0.6,
                                 stepSizeCoef = 4,
                                 delay = 10,
                                 assumeConvergence = 1000,
                                 trimSample = 40,
                                 maxiter = 3000))
      if(is.null(result)) next
      conditional <- result$conditional
      selected <- coordinates$selected[cluster$row]
      signal <- coordinates$signal[cluster$row]
      coordinatedat <- data.frame(conditional = conditional,
                                  observed = observed,
                                  signal = signal,
                                  selected = selected)
      coordinatedat$lCI[selected] <- result$coordinateCI[, 2]
      coordinatedat$uCI[selected] <- result$coordinateCI[, 1]

      try(profile <- optimizeSelected(observed, subCov, threshold,
                                      projected = mean(signal[selected]),
                                     stepRate = 0.6,
                                     stepSizeCoef = 4,
                                     delay = 10,
                                     assumeConvergence = 500,
                                     trimSample = 40,
                                     maxiter = 2000))

      naive <- mean(observed[selected])
      samp <- profile$sample
      profMeans <- rowMeans(samp[, selected, drop = FALSE])
      profPval <- 2 * min(mean(profMeans <= naive), mean(profMeans >= naive))

      e <- selected / sum(selected)
      naiveVar <- t(e) %*% subCov %*% e
      naiveSD <- sqrt(naiveVar)
      naiveCI <- naive + c(-1, 1) * 1.96 * naiveSD

      true <- mean(signal[selected])
      conditional <- mean(conditional[selected])
      lCI <- c(profPval, naiveCI[1], sort(result$meanCI)[1])
      uCI <- c(profPval, naiveCI[2],  sort(result$meanCI)[2])
      meanResult <- data.frame(type = c("true", "naive", "conditional"),
                               estimate = c(true, naive, conditional),
                               lCI = lCI, uCI = uCI)
      results[[m]][[3]] <- result
      results[[m]][[1]] <- coordinatedat
      results[[m]][[2]] <- meanResult
      results[[m]][[4]] <- profPval

      print(meanResult)

      iterCover <- iterCover + sum(selected) * (meanResult$lCI[3] < true & meanResult$uCI[3] > true)
      naiveCover <- naiveCover + sum(selected) * (meanResult$lCI[2] < true & meanResult$uCI[2] > true)
      profCover <- profCover + sum(selected) * (profPval > 0.05)
      weights <- weights + sum(selected)
    }

    simcover[rep, 1] <- sum(naiveCover) / sum(weights)
    simcover[rep, 2] <- sum(iterCover) / sum(weights)
    simcover[rep, 3] <- sum(profCover) / sum(weights)
    print(colMeans(simcover[1:rep, , drop = FALSE]))
    simresults[[rep]] <- results
  }

  return(simresults)
}

configurations <- expand.grid(snr = c(0.2, 0.8, 0.05),
                              rho = 0.8,
                              BHlevel = 0.1,
                              replications = 10)

set.seed(502)
system.time(simResults <- apply(configurations, 1, run.sim))
