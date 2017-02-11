library(selectiveMRI)

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

findLargestCluster <- function(x) {
  clustNum <- 1
  clustList <- list()
  cluster <- x[1]
  for(i in 2:(length(x))) {
    if(x[i] != x[i -1] + 1) {
      clustList[[clustNum]] <- cluster:x[i - 1]
      cluster <- x[i]
      clustNum <- clustNum + 1
    }
  }
  clustList[[clustNum]] <- cluster:x[length(x)]
  sizes <- sapply(clustList, length)
  return(clustList)
}

#set.seed(1231)
p <- 200
signalsig <- 0.05
signal <- 0.1
mu <- generateHMMmu(p, signal = signal, sig = signalsig)
mu <- predict(smooth.spline(mu))$y
rho <- 0.85
sigma <- rho^as.matrix(dist(1:p, diag = TRUE, upper = TRUE))

n <- 200
sample <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma)
cov <- var(sample)
testsd <- mean(diag(cov) / sqrt(n))
yhat <- colMeans(sample)
threshold <- 1.96 * testsd

plot(colMeans(sample), col = "red", type = "l")
lines(mu, type = "l")
abline(h = 0)
abline(h = c(-threshold, threshold), col = "grey")

select <- which(abs(yhat) > threshold)
clusters <- findLargestCluster(select)
means <- matrix(nrow = length(clusters), ncol = 4)
muEst <- yhat
for(i in 1:length(clusters)) {
  cluster <- clusters[[i]]
  if(length(cluster) <= 2) next
  print(round(i / length(clusters), 2))
  clusterPlus <- (min(cluster) - 1):(max(cluster) + 1)
  if(clusterPlus[1] == 0) clusterPlus <- clusterPlus[-1]
  if(clusterPlus[length(clusterPlus)] == p + 1) clusterPlus <- clusterPlus[-length(clusterPlus)]
  subCov <- cov[clusterPlus, clusterPlus] / (n)
  result <- optimizeSelected(yhat[clusterPlus], subCov, threshold, trimSample = 100,
                             maxiter = 3000)$estimates
  selected <- abs(yhat[clusterPlus]) > threshold
  conditional <- colMeans(result, floor(nrow(result) / 3):nrow(result))[selected]
  true <- mu[clusterPlus][selected]
  naive <- yhat[clusterPlus][selected]
  means[i, ] <- c(length(cluster), mean(true), mean(naive), mean(conditional))
  muEst[cluster] <- conditional
}
lines((muEst), col = "blue", type = "l")
means <- means[!is.na(means[, 1]), ]
print(means)


# Constrained -----------------
projectTo <- NULL
slack <- 10^-2
result <- optimizeSelected(yhat[clusterPlus], subCov, threshold,
                           projected = projectTo, quadraticSlack = slack,
                           lambdaStart = 50,
                           stepSizeCoef = 0.25,
                           stepRate = 0.5,
                           maxiter = 10^4^2,
                           assumeConvergence = 4000)$estimates
nsamp <- nrow(result)
conditional <- colMeans(result[floor(nsamp * 3 / 4):nsamp, selected])
mean(conditional)
c(mean(conditional^2), sum(selected) * projectTo^2 + slack)

require(dplyr)
require(reshape2)
require(ggplot2)
forplot <- melt(result[, selected])
names(forplot) <- c("iter", "variable", "mu")
ggplot(forplot) + geom_line(aes(x = iter, y = mu, col = factor(variable))) +
  theme_bw() + geom_hline(yintercept = 0)

