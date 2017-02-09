one.sided.pval <- function(x , threshold, Hnull = 0, sd = 1) {
  signx <- sign(x)
  x <- x * signx
  Hnull <- Hnull * signx
  pval <- 1 - ptruncnorm(x, a = abs(threshold), b = Inf, mean = Hnull)
  return(pval)
}

two.sided.probability <- function(x, threshold, Hnull = 0, sd = 1) {
  numerator <- numeric(length(x))
  numerator[x <= threshold] <- pnorm(x[x <= threshold], mean = Hnull, sd = sd)
  nLarger <- sum(x >= threshold)
  upper <- pnorm(-threshold, mean = Hnull, sd = sd) + pnorm(x[x >= threshold], mean = Hnull, sd = sd) -
    pnorm(threshold, mean = Hnull, sd = sd)
  numerator[x >= threshold] <- upper
  denominator <- 1 - pnorm(threshold, mean = Hnull, sd = sd) + pnorm(-threshold, mean = Hnull, sd = sd)
  return(numerator / denominator)
}

two.sided.pval <- function(x, threshold, Hnull = 0, sd = 1) {
  pval <- (1 - two.sided.probability(abs(x), abs(threshold), abs(Hnull), sd = sd))
  return(pval)
}

require(truncnorm)
threshold <- 1
reps <- 10^3
mu <- runif(reps, min = 0, max = 3.5)
ysign <- 1 - 2*rbinom(reps, 1, 1 - pnorm(threshold, mean = mu))
lthreshold <- rep(threshold, reps)
lthreshold[ysign == -1] <- -Inf
uthreshold <- rep(Inf, reps)
uthreshold[ysign == -1] <- -threshold
y <- rtruncnorm(reps, mean = mu, a = lthreshold, b = uthreshold)
onesided <- one.sided.pval(y, threshold)
twosided <- two.sided.pval(y, threshold)
plot(y, onesided, pch = ".", col = "blue")
points(y, twosided, col = "red", pch = ".")
abline(h = 0.05, col = "blue")
abline(h = 0.025, col = "red")
sum(twosided < 0.025)
sum(onesided < 0.05)


compute.onesided.mle <- function(y, threshold, sd = 1) {
  dens <- function(m) dtruncnorm(y, a = threshold, b = Inf, mean = m, sd = sd)
  fit <- optimize(dens, interval = c(-10, y), maximum = TRUE)
  return(fit$maximum)
}

mu <- seq(from = -5, to = 4, by = 0.02)
y <- qtruncnorm(0.5, a = 1.96, b = Inf, mean = mu, sd = 1)
plot(y, mu, type = "l", xlim = c(1.5, 4))
mle <- sapply(y, compute.onesided.mle, threshold)
lines(y, mle, col = "red")
abline(a = 0, b = 1, v = 1.96)
legend("bottomright", col = c("black", "red"), legend = c("mle", "median"), lty = 1)





