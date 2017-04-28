plotBRAIN <- function(coordinates, column, col = NULL) {
  for(ind in 1:K) {
    temp <- subset(coordinates, k == ind)
    signal <- temp[, column]
    imagemat <- matrix(nrow = I, ncol = J)
    for(l in 1:nrow(temp)) {
      imagemat[temp$i[l], temp$j[l]] <- signal[l]
    }

    if(is.null(col)) {
      (image(1:I, 1:J, imagemat, zlim = c(-max(abs(signal)), max(abs(signal))),
             main = ind))
    } else {
      (image(1:I, 1:J, imagemat, zlim = c(-max(abs(signal)), max(abs(signal))),
             main = ind, col = col))
    }
  }
}

findClusters <- function(coordinates) {
  columns <- which(names(coordinates) %in% c("i", "j", "k"))
  selected <- coordinates[coordinates$selected, ]
  graph <- as.matrix(dist(coordinates[, columns], method = "manhattan"))
  graph[graph > 1] <- 0
  clusterNumber <- 1
  clusters <- list()
  while(nrow(selected) > 0) {
    cluster <- selected[1, columns]
    toadd <- which(graph[selected$row[1], ] != 0)
    toVerify <- toadd[which(coordinates$selected[toadd])]
    cluster <- rbind(cluster, coordinates[setdiff(toadd, toVerify), columns])
    selected <- selected[-1, ]
    while(length(toVerify) > 0) {
      srow <- which(selected$row == toVerify[1])
      toVerify <- toVerify[-1]
      cluster <- rbind(cluster, selected[srow, columns])
      toadd <- which(graph[selected$row[srow], ] != 0)
      newVerify <- toadd[which(toadd %in% setdiff(selected$row, toVerify))]
      toVerify <- c(toVerify, newVerify)
      selected <- selected[-srow, ]
      cluster <- rbind(cluster, coordinates[setdiff(toadd, toVerify), columns])
    }

    cluster <- unique(cluster)
    cluster$row <- as.numeric(rownames(cluster))
    cluster$selected <- coordinates$selected[cluster$row]
    clusters[[clusterNumber]] <- cluster
    clusterNumber <- clusterNumber + 1
  }

  return(clusters)
}

# parameters + setup
# I <- 11
# J <- 10
# K <- 9
# rho <- 0.7
# coordinates <- expand.grid(i = 1:I, j = 1:J, k = 1:K)
# covariance <- rho^as.matrix(dist(coordinates[, 1:3], method = "euclidean",
#                                  diag = TRUE, upper = TRUE))
# covEigen <- eigen(covariance)
# sqrtCov <- covEigen$vectors %*% diag(sqrt(covEigen$values)) %*% t(covEigen$vectors)
# precision <- covEigen$vectors %*% diag((covEigen$values)^-1) %*% t(covEigen$vectors)
targetSnr <- 3.5

# Generating Signal ------------
coordinates <- expand.grid(i = 1:I, j = 1:J, k = 1:K)
signalProp <- 1
nnodes <- I * J * K
coordinates$row <- 1:nrow(coordinates)

mu <- sapply(c(I, J, K), function(x) rnorm(x))
mu <- apply(coordinates, 1, function(x) {
  sum(mu[[1]][1:x[1]]) + sum(mu[[2]][x[2]:J]) + sum(mu[[3]][1:x[3]])
})

mu <- mu - mean(mu)
coordinates$signal <- mu
par(mfrow = c(3, 3), mar = rep(2, 4))

location <- sapply(c(I, J, K), function(x) sample.int(x, 1))
s <- matrix(0.3, nrow = 3, ncol = 3)
diag(s) <- 1
s <- s*5
mu <- mvtnorm::dmvnorm(coordinates[, 1:3], mean = location, sigma = s)
mu <- mu * targetSnr / max(mu)
coordinates$signal <- mu

# Generating noise + data -----------------
noise <- rnorm(nnodes)
noise <- sqrtCov %*% noise
coordinates$noise <- noise
# snr <- var(coordinates$signal) / var(coordinates$noise)
# coordinates$signal <- coordinates$signal / sqrt(snr) * sqrt(targetSnr)
coordinates$observed <- coordinates$signal + coordinates$noise
par(mfrow = c(3, 3), mar = rep(2, 4))
plotBRAIN(coordinates, which(names(coordinates) == "signal"), col = rainbow(100))
plotBRAIN(coordinates, which(names(coordinates) == "noise"), col = rainbow(100))
plotBRAIN(coordinates, which(names(coordinates) == "observed"), col = rainbow(100))

# Univariate screening ----------------------
threshold <- 1.96
BHlevel <- 0.1
coordinates$zval <- coordinates$observed / sqrt(diag(covariance))
coordinates$pval <- 2 * pnorm(-abs(coordinates$zval))
coordinates$qval <- p.adjust(coordinates$pval, method = "BH")
par(mfrow = c(1, 1))
hist(coordinates$pval)
coordinates$selected <- coordinates$qval < BHlevel
# coordinates$selected <- abs(coordinates$observed) > threshold
par(mfrow = c(3, 3), mar = rep(2, 4))
plotBRAIN(coordinates, which(names(coordinates) == "signal"))
plotBRAIN(coordinates, which(names(coordinates) == "selected"))

# Inference  ----------------------------
clusters <- findClusters(coordinates)
sizes <- sapply(clusters, nrow)
cluster <- clusters[[1]]
cbind(coordinates$observed[cluster$row], cluster$selected)
threshold <- qnorm(BHlevel * sum(coordinates$selected) / nrow(coordinates) / 2,
                   lower.tail = FALSE)

results <- list()
pvals <- numeric(length(clusters))
par(mfrow = c(1, 1))
coordinates$estimate <- 0
coordinates$cluster <- 0
for(m in 1:length(clusters)) {
  results[[m]] <- list()
  cluster <- clusters[[m]]
  print(c(round(m / length(clusters), 2), nrow(cluster)))
  subCov <- covariance[cluster$row, cluster$row]
  observed <- coordinates$observed[cluster$row]
  selected <- coordinates$selected[cluster$row]
  signal <- coordinates$signal[cluster$row]
  try(result <- optimizeSelected(observed, subCov, threshold,
                                 projected = NULL,
                                 selected = selected,
                                 stepRate = 0.65,
                                 coordinates = cluster[, 1:3],
                                 tykohonovParam = NULL,
                                 tykohonovSlack = 1,
                                 stepSizeCoef = 1.5,
                                 delay = 20,
                                 assumeConvergence = 1800,
                                 trimSample = 50,
                                 maxiter = 2000,
                                 probMethod = "selected",
                                 init = observed,
                                 imputeBoundary = "neighbors"))
  #print(result$meanCI)
  samp <- rowMeans(result$sample[, selected, drop = FALSE])
  plot(density(samp), xlim = c(-5, 5))
  abline(v = mean(observed[selected]), col = "red")
  obsmean <- mean(observed[selected])
  pval <- 2 * min(mean(samp < obsmean), mean(samp > obsmean))
  pvals[m] <- pval
  print(c(pval = pval))

  cbind(colMeans(result$sample[, selected, drop = FALSE]), observed[selected])
  k <- 1
  plot(result$estimates[, selected, drop = FALSE][ ,k])
  abline(h = observed[selected][k])
  abline(h = signal[selected][k], col = "red")
  cbind(observed[selected], result$conditxional[selected], signal[selected])

  # try(truesamp <- optimizeSelected(observed, subCov, threshold,
  #                                 selected = selected,
  #                                 projected = mean(signal[selected]),
  #                                 stepRate = 0.6,
  #                                 coordinates = cluster[, 1:3],
  #                                 tykohonovParam = NULL,
  #                                 tykohonovSlack = 1,
  #                                 stepSizeCoef = 0,
  #                                 delay = 10,
  #                                 assumeConvergence = 2,
  #                                 trimSample = 50,
  #                                 maxiter = 1000,
  #                                 probMethod = "selected",
  #                                 init = observed,
  #                                 imputeBoundary = "neighbors"))
  # samp <- rowMeans(truesamp$sample[, selected, drop = FALSE])
  lines(density(samp), col = 'blue')
  abline(v = mean(signal[selected]), col = "green")
  abline(v = mean(rowMeans(result$sample[, selected, drop = FALSE])), col = "pink")
  abline(v = result$meanCI, col = "dark green")
  abline(v = mean(result$conditional[selected]), col = "orange")
  conditional <- result$conditional
  #print(mean(conditional[selected]))
  #print(mean(signal[selected]))
  selected <- coordinates$selected[cluster$row]
  signal <- coordinates$signal[cluster$row]
  coordinatedat <- data.frame(conditional = conditional,
                              observed = observed,
                              signal = signal,
                              selected = selected)
  coordinatedat$lCI[selected] <- result$coordinateCI[, 2]
  coordinatedat$uCI[selected] <- result$coordinateCI[, 1]
  results[[m]][[3]] <- result
  results[[m]][[1]] <- coordinatedat
  results[[m]][[2]] <- c(size = sum(selected),
                         conditional = mean(conditional[selected]),
                         observed = mean(observed[selected]),
                         signal = mean(signal[selected]),
                         lCI = sort(result$meanCI)[1], uCI = sort(result$meanCI)[2])
  print(results[[m]][[2]])
  coordinates$estimate[cluster$row[selected]] <- conditional[selected]
  coordinates$cluster[cluster$row[selected]] <- m
}

par(mfrow = c(3, 3))
plotBRAIN(coordinates, 5, col = rainbow(100))
plotBRAIN(coordinates, 12, col = rainbow(100))

coordinates[coordinates[, 12] != 0, c(13, 1:3, 7, 5, 12)]


