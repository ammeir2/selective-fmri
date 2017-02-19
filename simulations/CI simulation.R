generateHMMmu <- function(p, signal, sig = 0.1) {
  stateMat <- matrix(0.05, nrow = 3, ncol = 3)
  stateMat[1, ] <- c(0.95, 0.05, 0)
  stateMat[2, ] <- c(0.05, 0.9, 0.05)
  stateMat[3, ] <- c(0, 0.05, 0.95)
  #diag(stateMat) <- 0.9
  states <- rep(2, p)
  for(i in 2:p) {
    states[i] <- sample(1:3, 1, prob = stateMat[states[i - 1], ])
  }

  means <- c(-signal, 0 , signal)
  mu <- rnorm(p, mean = means[states], sig)
}

run.sim <- function(config) {
  n <- config[[1]]
  p <- config[[2]]
  rho <- config[[3]]
  signal <- config[[4]]

  sigma <- rho^as.matrix(dist(1:p, diag = TRUE, upper = TRUE))
  mu <- generateHMMmu(p, signal = signal, sig = signalsig)
  mu <- predict(smooth.spline(mu))$y

  n <- 200
  sample <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma)
  cov <- sigma #var(sample)
  testsd <- mean(diag(cov) / sqrt(n))
  yhat <- colMeans(sample)
  threshold <- 1.96 * testsd

  select <- which(abs(yhat) > threshold)
  clusters <- findLargestCluster(select)
  means <- matrix(nrow = length(clusters), ncol = 4)
  muEst <- yhat
  CIs <- list()
  for(i in 1:length(clusters)) {
    cluster <- clusters[[i]]
    if(length(cluster) <= 4) next
    cat(round(i / length(clusters), 2), " ")
    clusterPlus <- (min(cluster) - 1):(max(cluster) + 1)
    if(clusterPlus[1] == 0) clusterPlus <- clusterPlus[-1]
    if(clusterPlus[length(clusterPlus)] == p + 1) clusterPlus <- clusterPlus[-length(clusterPlus)]
    subCov <- cov[clusterPlus, clusterPlus] / (n)
    result <- optimizeSelected(yhat[clusterPlus], subCov, threshold,
                               stepRate = 0.6,
                               delay = 100,
                               assumeConvergence = 3000,
                               trimSample = 100,
                               maxiter = 7000)
    conditional <- result$conditional
    CIs[[i]] <- result$meanCI
    selected <- abs(yhat[clusterPlus]) > threshold
    true <- mu[clusterPlus][selected]
    naive <- yhat[clusterPlus][selected]
    conditional <- conditional[selected]
    means[i, ] <- c(length(cluster), mean(true), mean(naive), mean(conditional))
    muEst[cluster] <- conditional
  }
  print("")
  means <- means[!is.na(means[, 1]), ]
  CIs <- do.call("rbind", CIs)
  true <- means[, 2]
  print(cbind(CIs[, 1], true, CIs[, 2]))
  print(mean(sapply(1:nrow(CIs), function(i) CIs[i, 1] < true[i] & CIs[i, 2] > true[i])))

}






configurations <- expand.grid(n = c(50, 100, 200, 400), p = 200, rho = 0.85,
                              signal = 0.1, )
