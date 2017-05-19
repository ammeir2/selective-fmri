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
  colnames(simcover) <- c("cover", "power")
  for(rep in 1:replications) {
    selected <- FALSE
    maxsize <- 0
    while(is.na(maxsize) | maxsize < 2) {
      coordinates <- generateArrayData3D(dims, sqrtCov, snr, spread)
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

    print(c(rep = rep, rho = rho, BHlevel = BHlevel, snr = snr, spread = spread))
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
        } else if (methodind == 3) {
          method <- "naive"
          results[[slot]] <- list()
          results[[slot]][[1]] <- c(snr = snr, spread = spread, method = methodind,
                                    rho = rho, BHlevel = BHlevel,
                                    size = sum(selected), true = true)
          e <- rep(1 / sum(selected), sum(selected))
          naiveVar <- as.numeric(t(e) %*% subCov[selected, selected] %*% e)
          naiveSD <- sqrt(naiveVar)
          zval <- naive / naiveSD
          pvalue <- pnorm(-abs(zval))
          profPval <- 2 * min(pnorm(naive, signal, naiveSD), pnorm(naive, signal, naiveSD, lower.tail = FALSE))
          profResult <- c(true = mean(signal[selected]), profPval = profPval, pvalue = pvalue)
          results[[slot]][[2]] <- profResult
          results[[slot]][[3]] <- c(conditional = NA, naive = naive, true = true)
          slot <- slot + 1
          next
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
                              BHlevel = c(0.1),
                              replications = 30)

system.time(simResults <- apply(configurations, 1, run.sim))
save(simResults, file = "simulations/results/May 12J naive.Robj")

load(file = "simulations/results/May 12 naive.Robj")
simResults12 <- simResults
load(file = "simulations/results/May 12B naive.Robj")
simResults12B <- simResults
load(file = "simulations/results/May 12C naive.Robj")
simResults12C <- simResults
load(file = "simulations/results/May 12D naive.Robj")
simResults12D <- simResults
load(file = "simulations/results/May 12E naive.Robj")
simResults12E <- simResults
load(file = "simulations/results/May 12F naive.Robj")
simResults12F <- simResults
load(file = "simulations/results/May 12G naive.Robj")
simResults12G <- simResults
load(file = "simulations/results/May 12H naive.Robj")
simResults12H <- simResults
load(file = "simulations/results/May 12I naive.Robj")
simResults12I <- simResults
load(file = "simulations/results/May 12J naive.Robj")
simResults12J <- simResults



simResults <- c(simResults12, simResults12B, simResults12C,
                simResults12D, simResults12E, simResults12F,
                simResults12G, simResults12H, simResults12I,
                simResults12J)

# Processing ------------------------
library(ggplot2)
library(dplyr)
library(reshape2)
resultList <- list()
len <- 1
for(i in 1:length(simResults)) {
  for(j in 1:length(simResults[[i]])) {
    nonEmpty <- sapply(simResults[[i]][[j]], length) > 0
    iterResult <- data.frame(matrix(nrow = sum(nonEmpty), ncol = 13))
    iterList <- simResults[[i]][[j]]
    row <- 1
    exp <- j * runif(1)
    for(k in which(nonEmpty)) {
      method <- iterList[[k]][[1]][[3]]
      if(method == 1) {
        method <- "onesided"
      } else if(method == 2) {
        method <- "selected"
      } else if (method == 3) {
        method <- "naive"
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
                             cover = as.numeric(iterList[[k]][[2]][[2]]),
                             naive = iterList[[k]][[3]][[2]],
                             cond = iterList[[k]][[3]][[1]])
      row <- row + 1
    }
    iterResult <- data.frame(iterResult)
    names(iterResult) <- c("experiment", "cluster", "snr", "spread", "method", "rho",
                           "pthreshold", "size", "true", "pval", "cover", "naive", "cond")
    resultList[[len]] <- iterResult
    len <- len + 1
  }
}
results <- do.call("rbind", resultList)

results$cover <- as.numeric(as.character(results$cover))
results$pval <- as.numeric(as.character(results$pval))
results$true <- as.numeric(as.character(results$true))
results$snr <- as.numeric(as.character(results$snr))
results$experiment <- as.numeric(as.character(results$experiment))
results$size <- as.numeric(as.character(results$size))
results$cond <- as.numeric(as.character(results$cond))
results$naive <- as.numeric(as.character(results$naive))


# Plot by snr --------------------
level <- 0.05
cover <- summarize(group_by(results, experiment, snr, spread, method, rho, pthreshold),
                   cover = sum((cover > level) * size, na.rm = TRUE) / sum(size[!is.na(cover > level)]),
                   power = pval[which.max(true)] < level)

cover <- summarize(group_by(cover, snr, spread, method, rho, pthreshold),
                   coversd = sd(cover) / sqrt(length(cover)),
                   cover = mean(cover),
                   powersd = sd(power, na.rm = TRUE) / sqrt(length(power)),
                   power = mean(power, na.rm = TRUE))


quantile <- qnorm(1 - 0.05 / nrow(cover))
ciwidth <- 0.1
coverplot <- ggplot(cover) +
  geom_line(aes(x = snr, y = cover, col = method, linetype = method)) +
  geom_point(aes(x = snr, y = cover, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_hline(yintercept = 1 - level) + theme_bw() +
  geom_segment(aes(y = pmax(cover - coversd * quantile, 0), yend = pmin(cover + coversd * quantile, 1),
                   x = snr, xend = snr, col = method, linetype = method)) +
  geom_segment(aes(y = pmax(cover - coversd * quantile, 0), yend = pmax(cover - coversd * quantile, 0),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  geom_segment(aes(y = pmin(cover + coversd * quantile, 1), yend = pmin(cover + coversd * quantile, 1),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  ggtitle("Coverage Rate") + theme(legend.position = "none") +
  scale_colour_brewer(palette = "Set1")
coverplot

powerplot <- ggplot(cover) +
  geom_line(aes(x = snr, y = power, col = method, linetype = method)) +
  geom_point(aes(x = snr, y = power, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_hline(yintercept = level) + theme_bw() +
  geom_segment(aes(y = pmax(power - powersd * quantile, 0), yend = pmin(power + powersd * quantile, 1),
                   x = snr, xend = snr, col = method, linetype = method)) +
  geom_segment(aes(y = pmax(power - powersd * quantile, 0), yend = pmax(power - powersd * quantile, 0),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  geom_segment(aes(y = pmin(power + powersd * quantile, 1), yend = pmin(power + powersd * quantile, 1),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  ylim(0, 1) + xlim(1 - ciwidth, 4 + ciwidth) +
  ggtitle("Power") + theme(legend.position = "none")
powerplot

# pdf("figures/ciplots.pdf",pagecentre=T, width=8, height=2.5 ,paper = "special")
# gridExtra::grid.arrange(coverplot, powerplot, nrow = 1)
# dev.off()

# Plot by true signal ---------------------------
# temp <- results
# ncut <- 4
# temp$true <- cut(temp$true, ncut, labels = FALSE)
# map <- cbind(1:ncut, seq(from = 0, to = max(results$true), length.out = ncut))
# temp$true <- map[temp$true, 2]
# cover <- summarize(group_by(temp, experiment, true, spread, method, rho, BHlevel),
#                    cover = sum((cover > level) * size, na.rm = TRUE) / sum(size[!is.na(cover > level)]),
#                    power = pval[which.max(true)] < level)
# cover <- summarize(group_by(cover, true, spread, method, rho, BHlevel),
#                    coversd = sd(cover) / sqrt(length(cover)),
#                    cover = mean(cover),
#                    powersd = sd(power, na.rm = TRUE) / sqrt(length(power)),
#                    power = mean(power, na.rm = TRUE))
#
#
# quantile <- qnorm(1 - 0.05 / nrow(cover))
# cover$true <- as.numeric(cover$true)
# ggplot(cover) +
#   geom_line(aes(x = true, y = cover, col = method)) +
#   facet_grid(rho ~ BHlevel, labeller = "label_both") +
#   geom_hline(yintercept = 1 - level) + theme_bw() +
#   geom_segment(aes(y = pmax(cover - coversd * quantile, 0), yend = pmin(cover + coversd * quantile, 1),
#                    x = true, xend = true, col = method), linetype = 2) +
#   ggtitle("Coverage Rate for Profile CIs")
#
# ggplot(cover) +
#   geom_line(aes(x = true, y = power, col = method)) +
#   facet_grid(rho ~ BHlevel, labeller = "label_both") +
#   geom_hline(yintercept = level) + theme_bw() +
#   geom_segment(aes(y = pmax(power - powersd * quantile, 0), yend = pmin(power + powersd * quantile, 1),
#                    x = true, xend = true, col = method), linetype = 2) +
#   ylim(0, 1) +
#   ggtitle("Power for Profile CIs")

# Estimation Error ---------------------
naive <- subset(results, method == "selected")
naive <- summarize(group_by(naive, snr, pthreshold, rho, experiment),
                  mse = weighted.mean(sqrt((naive - true)^2), size),
                  bias = weighted.mean(sign(naive) * (naive - true), size))
naive <- summarize(group_by(naive, snr, pthreshold, rho),
                   msesd = sd(mse) / sqrt(length(mse)), mse = mean(mse),
                   biassd = sd(bias) / sqrt(length(bias)), bias = mean(bias))
naive$method <- "naive"

cond <- subset(results, method == "selected")
cond <- summarize(group_by(cond, snr, pthreshold, rho, experiment),
                   mse = weighted.mean(sqrt((cond - true)^2), size),
                   bias = weighted.mean(sign(naive)*(cond - true), size))
cond <- summarize(group_by(cond, snr, pthreshold, rho),
                  msesd = sd(mse) / sqrt(length(mse)), mse = mean(mse),
                  biassd = sd(bias) / sqrt(length(bias)), bias = mean(bias))
cond$method <- "conditional"

forplot <- rbind(naive, cond)
mse <- ggplot(forplot, aes(x = snr, y = mse, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_line(aes(linetype = method)) + geom_point() +
  geom_segment(aes(y = mse - 2 * msesd, yend = mse + 2 * msesd,
                   x = snr, xend = snr, col = method), linetype = 2) +
  theme_bw() + ggtitle("Root Mean Squared Error") +
  theme(legend.position = "none")
mse

bias <- ggplot(forplot, aes(x = snr, y = bias, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_line(aes(linetype = method)) + geom_point() +
  geom_segment(aes(y = bias - 2 * biassd, yend = bias + 2 * biassd,
                   x = snr, xend = snr, col = method), linetype = 2) +
  geom_hline(yintercept = 0) + theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Bias")
bias

# pdf("figures/estplots.pdf",pagecentre=T, width=8, height=2.5 ,paper = "special")
# gridExtra::grid.arrange(mse, bias, nrow = 1)
# dev.off()




