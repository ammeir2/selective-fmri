#library(selective-fmri)
run.sim <- function(config) {
  snr <- config[[1]]
  spread <- config[[2]]
  rho <- config[[3]]
  BHlevel <- config[[4]]
  replications <- config[[5]]

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
  simcover <- matrix(nrow = replications, ncol = 2)
  colnames(simcover) <- c("prof", "better")
  for(rep in 1:replications) {
    selected <- FALSE
    maxsize <- 0
    while(is.na(maxsize) | maxsize < 2) {
      coordinates <- generateArrayData3D(dims, sqrtCov, snr, spread)
      coordinates$zval <- coordinates$observed / sqrt(diag(covariance))
      coordinates$pval <- 2 * pnorm(-abs(coordinates$zval))
      coordinates$qval <- p.adjust(coordinates$pval, method = "BH")
      coordinates$selected <- coordinates$qval < BHlevel
      selected <- coordinates$selected
      if(sum(selected) >= 2) {
        clusters <- findClusters(coordinates)
        maxsize <- max(sapply(clusters, function(x) sum(x$selected)))
      }
    }

    print(c(rep = rep, rho = rho, BHlevel = BHlevel, snr = snr, spread = spread))
    print(c(nselected = sum(selected)))
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
      selected <- coordinates$selected[cluster$row]
      observed <- coordinates$observed[cluster$row]
      signal <- coordinates$signal[cluster$row]

      #print(rbind(observed[!selected], signal[!selected]))

      if(sum(cluster$selected) == 1) next
      print(c(round(m / length(clusters), 2), nrow(cluster)))
      subCov <- covariance[cluster$row, cluster$row, drop = FALSE]

      try(result <- optimizeSelected(observed, subCov, threshold,
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
                                      maxiter = 1005,
                                     probMethod = "onesided",
                                      init = rep(0, length(observed)),
                                     imputeBoundary = "neighbors"))
      if(is.null(result)) next
      naive <- mean(observed[selected])
      samp <- result$sample
      profMeans <- NULL
      try(profMeans <- rowMeans(samp[, selected, drop = FALSE]))
      if(is.null(profMeans)) next
      betterPval <-  min(mean(naive < profMeans), mean(naive > profMeans))
      # if(naive > 0) {
      #   betterPval <- mean(naive < profMeans)
      # } else {
      #   betterPval <- mean(naive > profMeans)
      # }

      obsratio <- abs(mean(observed[selected]) / mean(signal[selected]))
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
                                     maxiter = 1800,
                                     probMethod = "onesided",
                                     init = observed,
                                     imputeBoundary = "neighbors"))

      samp <- profile$sample
      profMeans <- NULL
      try(profMeans <- rowMeans(samp[, selected, drop = FALSE]))
      if(is.null(profMeans)) next
      profPval <-  min(mean(naive < profMeans), mean(naive > profMeans))
      # if(naive > 0) {
      #   profPval <- mean(naive < profMeans)
      # } else {
      #   profPval <- mean(naive > profMeans)
      # }

      true <- mean(signal[selected])
      profResult <- c(true = mean(signal[selected]), profPVAL = profPval, betterPVAL = betterPval)
      results[[m]][[1]] <- c(snr = snr, spread = spread, rho = rho, BHlevel = BHlevel, size = sum(selected),
                             true = true)
      results[[m]][[2]] <- profResult

      print(profResult)

      weight <- 1
      iterCover <- iterCover + weight * (betterPval < 0.05)
      profCover <- profCover + weight * (profPval > 0.05)
      weights <- weights + weight

      plot(density((rowMeans(profile$sample[, selected]))), xlim = c(-4.5, 4.5))
      abline(v = mean(profile$conditional[selected]))
      abline(v = mean(result$conditional[selected]), col = "blue")
      abline(v = naive, col = "red")
      lines(density((rowMeans(result$sample[, selected]))), col = "blue")
      abline(v = c(threshold, - threshold), col = "grey")
      # print(c(threshold, min(abs(profile$sample[-(1:10), selected]))))
      # print(c(threshold, min(abs(result$sample[-(1:10), selected]))))
      # fff <- 1
    }

    simcover[rep, 1] <- sum(profCover) / sum(weights)
    simcover[rep, 2] <- sum(iterCover) / sum(weights)
    print(c(iterprof = profCover, itermle = iterCover) / weights)
    print(colMeans(simcover[1:rep, , drop = FALSE]))
    simresults[[rep]] <- results
  }

  return(simresults)
}

configurations <- expand.grid(snr = c(3, 2.5, 2, 1.5, 1),
                              spread = c(2, 1),
                              rho = c(0.5, 0.75),
                              BHlevel = 0.1,
                              replications = 50)

set.seed(510)
system.time(simResults <- apply(configurations, 1, run.sim))
#save(simResults, file = "simulations/results/Apr 18 twosided w power.Robj")
#load(file = "simulations/results/Apr 17 twosided w power.Robj")

# Processing ------------------------
library(ggplot2)
library(dplyr)
library(reshape2)
resultList <- list()
len <- 1
for(i in 1:length(simResults)) {
  for(j in 1:length(simResults[[i]])) {
    nonEmpty <- sapply(simResults[[i]][[j]], length) > 0
    iterResult <- matrix(nrow = sum(nonEmpty), ncol = 10)
    iterList <- simResults[[i]][[j]]
    row <- 1
    for(k in which(nonEmpty)) {
      iterResult[row, ] <- c(experiment = j,
                             cluster = row,
                             snr = iterList[[k]][[1]][[1]],
                             spread = iterList[[k]][[1]][[2]],
                             rho = iterList[[k]][[1]][[3]],
                             BHlevel = iterList[[k]][[1]][[4]],
                             size = iterList[[k]][[1]][[5]],
                             true = iterList[[k]][[1]][[6]],
                             pval = iterList[[k]][[2]][[3]],
                             cover = iterList[[k]][[2]][[2]])
      row <- row + 1
    }
    iterResult <- data.frame(iterResult)
    names(iterResult) <- c("experiment", "cluster", "snr", "spread", "rho",
                           "BHlevel", "size", "true", "pval", "cover")
    resultList[[len]] <- iterResult
    len <- len + 1
  }
}
results <- do.call("rbind", resultList)

cover <- summarize(group_by(results, experiment, snr, spread, rho),
                   cover = mean(cover > 0.05),
                   power = mean(pval[true > 0.05] < 0.05, na.rm = TRUE))
cover <- summarize(group_by(cover, snr, spread, rho),
                   coversd = sd(cover) / sqrt(length(cover)),
                   cover = mean(cover),
                   powersd = sd(power, na.rm = TRUE) / sqrt(length(power)),
                   power = mean(power, na.rm = TRUE))
