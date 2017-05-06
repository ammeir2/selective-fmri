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
    slot <- 1
    for(m in 1:length(clusters)) {
      cluster <- clusters[[m]]
      cluster <- subset(cluster, !is.na(cluster$selected))
      selected <- coordinates$selected[cluster$row]
      observed <- coordinates$observed[cluster$row]
      signal <- coordinates$signal[cluster$row]
      for(methodind in 1:2) {
        if(methodind == 1) {
          method <- "onesided"
        } else if(methodind == 2) {
          method <- "selected"
        }
        results[[slot]] <- list()

        #print(rbind(observed[!selected], signal[!selected]))
        #if(sum(cluster$selected) == 1) next
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
                                       probMethod = method,
                                       init = rep(0, length(observed)),
                                       imputeBoundary = "neighbors"))
        if(is.null(result)) next
        naive <- mean(observed[selected])
        samp <- result$sample
        profMeans <- NULL
        try(profMeans <- rowMeans(samp[, selected, drop = FALSE]))
        if(is.null(profMeans)) next
        betterPval <-  min(mean(naive < profMeans), mean(naive > profMeans))

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
                                        probMethod = method,
                                        init = observed,
                                        imputeBoundary = "neighbors"))

        samp <- profile$sample
        profMeans <- NULL
        try(profMeans <- rowMeans(samp[, selected, drop = FALSE]))
        if(is.null(profMeans)) next
        profPval <-  min(mean(naive < profMeans), mean(naive > profMeans))

        true <- mean(signal[selected])
        profResult <- c(true = mean(signal[selected]), profPVAL = profPval, betterPVAL = betterPval)
        results[[slot]][[1]] <- c(snr = snr, spread = spread, method = methodind,
                               rho = rho, BHlevel = BHlevel,
                               size = sum(selected), true = true)
        results[[slot]][[2]] <- profResult

        print(profResult)

        weight <- 1
        iterCover <- iterCover + weight * (betterPval < 0.05)
        profCover <- profCover + weight * (profPval > 0.05)
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

    simcover[rep, 1] <- sum(profCover) / sum(weights)
    simcover[rep, 2] <- sum(iterCover) / sum(weights)
    print(c(iterprof = profCover, itermle = iterCover) / weights)
    print(colMeans(simcover[1:rep, , drop = FALSE]))
    simresults[[rep]] <- results
  }

  return(simresults)
}

configurations <- expand.grid(snr = c(3.5, 3, 2.5, 2, 1.5, 1),
                              spread = c(2, 1),
                              rho = c(0.5, 0.75),
                              BHlevel = 0.1,
                              replications = 50)

system.time(simResults <- apply(configurations, 1, run.sim))
#save(simResults, file = "simulations/results/May 5C both.Robj")

load(file = "simulations/results/May 4 both.Robj")
simResults4 <- simResults
load(file = "simulations/results/May 4B both.Robj")
simResults4B <- simResults
load(file = "simulations/results/May 4C both.Robj")
simResults4C <- simResults
load(file = "simulations/results/May 5 both.Robj")
simResults5 <- simResults
load(file = "simulations/results/May 5B both.Robj")
simResults5B <- simResults
load(file = "simulations/results/May 5B both.Robj")
simResults5C <- simResults


simResults <- c(simResults4, simResults4B, simResults4C,
                simResults5, simResults5B, simResults5C)

# Processing ------------------------
library(ggplot2)
library(dplyr)
library(reshape2)
resultList <- list()
len <- 1
for(i in 1:length(simResults)) {
  for(j in 1:length(simResults[[i]])) {
    nonEmpty <- sapply(simResults[[i]][[j]], length) > 0
    iterResult <- data.frame(matrix(nrow = sum(nonEmpty), ncol = 11))
    iterList <- simResults[[i]][[j]]
    row <- 1
    exp <- j * runif(1)
    for(k in which(nonEmpty)) {
      method <- iterList[[k]][[1]][[3]]
      if(method == 1) {
        method <- "onesided"
      } else if(method == 2) {
        method <- "selected"
      }
      iterResult[row, ] <- c(experiment = exp,
                             cluster = row,
                             snr = as.numeric(iterList[[k]][[1]][[1]]),
                             spread = as.numeric(iterList[[k]][[1]][[2]]),
                             method = method,
                             rho = as.numeric(iterList[[k]][[1]][[4]]),
                             BHlevel = as.numeric(iterList[[k]][[1]][[5]]),
                             size = as.numeric(iterList[[k]][[1]][[6]]),
                             true = as.numeric(iterList[[k]][[1]][[7]]),
                             pval = as.numeric(iterList[[k]][[2]][[3]]),
                             cover = as.numeric(iterList[[k]][[2]][[2]]))
      row <- row + 1
    }
    iterResult <- data.frame(iterResult)
    names(iterResult) <- c("experiment", "cluster", "snr", "spread", "method", "rho",
                           "BHlevel", "size", "true", "pval", "cover")
    resultList[[len]] <- iterResult
    len <- len + 1
  }
}
results <- do.call("rbind", resultList)
two <- readRDS(file = "simulations/results/may4two.Robj")
one <- readRDS(file = "simulations/results/may4one.Robj")
results <- rbind(one, two, results)
results$cover <- as.numeric(as.character(results$cover))
results$pval <- as.numeric(as.character(results$pval))
results$true <- as.numeric(as.character(results$true))
results$snr <- as.numeric(as.character(results$snr))
results$experiment <- as.numeric(as.character(results$experiment))
results$size <- as.numeric(as.character(results$size))

level <- 0.1 / 2
cover <- summarize(group_by(results, experiment, snr, spread, method, rho),
                   cover = sum((cover > level) * size, na.rm = TRUE) / sum(size[!is.na(cover > level)]),
                   power = mean(pval[true > 0.05] < level, w = size, na.rm = TRUE))
cover <- summarize(group_by(cover, snr, spread, method, rho),
                   coversd = sd(cover) / sqrt(length(cover)),
                   cover = mean(cover),
                   powersd = sd(power, na.rm = TRUE) / sqrt(length(power)),
                   power = mean(power, na.rm = TRUE))


quantile <- qnorm(1 - 0.05 / nrow(cover))
ggplot(cover) +
  geom_line(aes(x = snr, y = cover, col = method)) +
  facet_grid(rho ~ spread, labeller = "label_both") +
  geom_hline(yintercept = 1 - level * 2) + theme_bw() +
  geom_segment(aes(y = pmax(cover - coversd * quantile, 0), yend = pmin(cover + coversd * quantile, 1),
                   x = snr, xend = snr, col = method), linetype = 2) +
  ggtitle("Coverage Rate for Profile CIs")

ggplot(cover) +
  geom_line(aes(x = snr, y = power, col = method)) +
  facet_grid(rho ~ spread, labeller = "label_both") +
  geom_hline(yintercept = level) + theme_bw() +
  geom_segment(aes(y = pmax(power - powersd * quantile, 0), yend = pmin(power + powersd * quantile, 1),
                   x = snr, xend = snr, col = method), linetype = 2) +
  ylim(0, 1) +
  ggtitle("Power for Profile CIs")

