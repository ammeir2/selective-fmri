pseq <- seq(from = 1, to = 200, by = 10)
fullMu <- rnorm(pmax)
rho <- 0.3
benchmark <- numeric(length(pseq))
reps <- 10
for(j in 1:length(pseq)) {
  p <- pseq[j]
  times <- numeric(reps)
  for(i in 1:reps) {
    rho <- runif(1, min = 0.1, max = 0.9)
    a <- rexp(p)
    b <- rep(Inf, p)
    mu <- rnorm(p)
    sigma <- matrix(rho, nrow = p, ncol = p)
    diag(sigma) <- 1
    times[i] <- system.time(mvtnorm::pmvnorm(lower = a, upper = b, mean = mu, sigma = sigma))[3]
  }
  benchmark[j] <- mean(times)
  print(c(p, benchmark[j]))
}

plot(pseq, benchmark, type = "b")
abline(a = 0, b = 1)
lm(benchmark ~ pseq)
