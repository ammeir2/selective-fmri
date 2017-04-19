#library(selective-fmri)
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
  simcover <- matrix(nrow = replications, ncol = 2)
  colnames(simcover) <- c("prof", "better")
  for(rep in 1:replications) {
    selected <- FALSE
    maxsize <- 0
    while(is.na(maxsize) | maxsize < 2) {
      coordinates <- generateArrayData3D(dims, sqrtCov, snr)
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

    print(c(rep = rep, rho = rho, BHlevel = BHlevel, snr = snr))
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

      result <- NULL
      try(result <- optimizeSelected(observed, subCov, threshold,
                                      projected = 0,
                                     selected = selected,
                                      stepRate = 0.6,
                                      coordinates = cluster[, 1:3],
                                     tykohonovParam = NULL,
                                     tykohonovSlack = 0.00001,
                                      stepSizeCoef = 0,
                                      delay = 1,
                                      assumeConvergence = 2,
                                      trimSample = 50,
                                      maxiter = 2500,
                                     probMethod = "all",
                                      init = rep(0, length(observed)),
                                     imputeBoundary = "neighbors"))
      if(is.null(result)) next
      naive <- mean(observed[selected])
      samp <- result$sample
      profMeans <- rowMeans(samp[, selected, drop = FALSE])
      betterPval <- 2 * min(mean(naive < profMeans), mean(naive > profMeans))
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
                                     tykohonovSlack = min(1, obsratio),
                                     stepSizeCoef = 4,
                                     delay = 10,
                                     assumeConvergence = 500,
                                     trimSample = 50,
                                     maxiter = 2500,
                                     probMethod = "all",
                                     init = observed,
                                     imputeBoundary = "neighbors"))

      samp <- profile$sample
      profMeans <- rowMeans(samp[, selected, drop = FALSE])
      profPval <- 2 * min(mean(naive < profMeans), mean(naive > profMeans))
      # if(naive > 0) {
      #   profPval <- mean(naive < profMeans)
      # } else {
      #   profPval <- mean(naive > profMeans)
      # }

      true <- mean(signal[selected])
      profResult <- c(true = mean(signal[selected]), profPVAL = profPval, betterPVAL = betterPval)
      results[[m]][[1]] <- c(snr = snr, rho = rho, BHlevel = BHlevel, size = sum(selected))
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
    }

    simcover[rep, 1] <- sum(profCover) / sum(weights)
    simcover[rep, 2] <- sum(iterCover) / sum(weights)
    print(c(iterprof = profCover, itermle = iterCover) / weights)
    print(colMeans(simcover[1:rep, , drop = FALSE]))
    simresults[[rep]] <- results
  }

  return(simresults)
}

configurations <- expand.grid(snr = c(0.5, 0.25, 0.1),
                              rho = c(0.5, 0.75),
                              BHlevel = 0.1,
                              replications = 50)

#set.seed(510)
system.time(simResults <- apply(configurations, 1, run.sim))
save(simResults, file = "simulations/results/Apr 18 twosided w power.Robj")
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
    iterResult <- matrix(nrow = sum(nonEmpty), ncol = 8)
    iterList <- simResults[[i]][[j]]
    row <- 1
    for(k in which(nonEmpty)) {
      iterResult[row, ] <- c(experiment = j,
                             cluster = row,
                             snr = iterList[[k]][[1]][[1]],
                             rho = iterList[[k]][[1]][[2]],
                             BHlevel = iterList[[k]][[1]][[3]],
                             size = iterList[[k]][[1]][[4]],
                             pval = iterList[[k]][[2]][[3]],
                             cover = iterList[[k]][[2]][[2]])
      row <- row + 1
    }
    iterResult <- data.frame(iterResult)
    names(iterResult) <- c("experiment", "cluster", "snr", "rho",
                           "BHlevel", "size", "pval", "cover")
    resultList[[len]] <- iterResult
    len <- len + 1
  }
}
results <- do.call("rbind", resultList)

cover <- summarize(group_by(results, experiment, snr, rho),
                   cover = mean(cover > 0.05),
                   power = mean(pval < 0.05))
cover <- summarize(group_by(cover, snr, rho),
                   coversd = sd(cover) / sqrt(length(cover)),
                   cover = mean(cover),
                   powersd = sd(power) / sqrt(length(power)),
                   power = mean(power))
