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

      unselected <- which(!selected)
      distances <- as.matrix(dist(cluster[, 1:3]))
      diag(distances) <- Inf
      neighbors <- cbind(unselected, apply(distances[unselected, ], 1, function(x) which(selected)[which.min(x[selected])[1]]))

      result <- NULL
      try(result <- optimizeSelected(observed, subCov, threshold,
                                      projected = 0,
                                      stepRate = 0.6,
                                     quadraticSlack = 3,
                                      stepSizeCoef = 4,
                                      delay = 1,
                                      assumeConvergence = 500,
                                      trimSample = 50,
                                      maxiter = 2500,
                                     probMethod = "onesided",
                                      init = observed,
                                     imputeBoundary = TRUE,
                                     neighbors = neighbors))
      if(is.null(result)) next
      naive <- mean(observed[selected])
      samp <- result$sample
      profMeans <- rowMeans(samp[, selected, drop = FALSE])
      betterPval <- 2 * min(mean(naive < profMeans), mean(naive > profMeans))
      if(naive > 0) {
        betterPval <- mean(naive < profMeans)
      } else {
        betterPval <- mean(naive > profMeans)
      }


      try(profile <- optimizeSelected(observed, subCov, threshold,
                                      projected = mean(signal[selected]),
                                     stepRate = 0.6,
                                     quadraticSlack = 3,
                                     stepSizeCoef = 4,
                                     delay = 1,
                                     assumeConvergence = 500,
                                     trimSample = 50,
                                     maxiter = 2500,
                                     probMethod = "onesided",
                                     init = observed,
                                     imputeBoundary = TRUE,
                                     neighbors = neighbors))

      samp <- profile$sample
      profMeans <- rowMeans(samp[, selected, drop = FALSE])
      profPval <- 2 * min(mean(naive < profMeans), mean(naive > profMeans))
      if(naive > 0) {
        profPval <- mean(naive < profMeans)
      } else {
        profPval <- mean(naive > profMeans)
      }

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

configurations <- expand.grid(snr = c(2, 0.5, 0.25),
                              rho = c(0.5, 0.75),
                              BHlevel = 0.1,
                              replications = 50)

set.seed(502)
system.time(simResults <- apply(configurations, 1, run.sim))
save(simResults, file = "simulations/results/Mar 27 oneside w power.Robj")

# One sided  w oracle --- load(file = "simulations/results/Mar 27.Robj")
# Two Sided  w oracle --- load(file = "simulations/results/Mar 26D.Robj")

# Processing ------------------------
library(ggplot2)
library(dplyr)
library(reshape2)
coverList <- list()
len <- 1
for(i in 1:length(simResults)) {
  setting <- simResults[[i]]
  for(j in 1:length(setting)) {
    iter <- setting[[j]]
    for(k in 1:length(iter)) {
      sim <- iter[[k]]
      if(length(sim) == 0) next
      params <- sim[[1]]
      pvals <- sim[[2]]
      coverList[[len]] <- c(iter = j, params, pvals)
      len <- len + 1
    }
  }
}
cover <- data.frame(do.call("rbind", coverList))
cover <- melt(cover, id.vars = c("iter", "snr", "rho", "BHlevel", "size"))
names(cover)[6:7] <- c("method", "pval")
cover <- summarize(group_by(cover, snr, rho, BHlevel, method, iter),
                   weighted = sum(size * (pval > 0.05)) / sum(size),
                   unweighted = mean(pval > 0.05))
cover <- melt(cover, id.vars = c("snr", "rho", "BHlevel", "method", "iter"))
names(cover)[6:7] <- c("weight", "cover")
meancover <- summarize(group_by(cover, snr, rho, BHlevel, method, weight),
                       cover = mean(cover))
ggplot(cover) +
  geom_boxplot(aes(x = method, y = cover, col = weight)) +
  facet_grid(snr ~ rho, labeller = "label_both") + theme_bw() +
  geom_hline(yintercept = 0.95) +
  geom_point(data = meancover, aes(x = method, y = cover, col = weight),
             size = 3, alpha = 0.5)

subset(data.frame(meancover), weight == "unweighted")
