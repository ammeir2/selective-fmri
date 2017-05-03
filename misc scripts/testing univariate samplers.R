# First example --------------
mu <- 0.2
lower <- -Inf
upper <- -1
sd <- 0.5
reps <- 10^4
extreme <- numeric(reps)
standard <- numeric(reps)
benchmark <- truncnorm::rtruncnorm(reps, lower, upper, mu, sd)
for(i in 1:reps) {
  extreme[i] <- sampleExtreme(mu, sd, lower, upper)
  standard[i] <- sampleUnivTruncNorm(mu, sd, lower, upper)
}
plot(density(benchmark))
lines(density(standard), col = "blue")
lines(density(extreme), col = "red")
