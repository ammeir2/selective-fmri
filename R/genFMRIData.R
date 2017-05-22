residualData3D <- function(dat, grp_size, targetSnr ,spread = 1) {

  dims = dim(dat)

  # Generating signal
  coordinates <- expand.grid(i=1:dims[1],j=1:dims[2],k=1:dims[3])
  nnodes <- prod(dims[1:3])
  coordinates$row = 1:nrow(coordinates)

  location <- sapply(dims[1:3], function(x) sample.int(x,1))
  s <- matrix(0.3, nrow = 3, ncol = 3)
  diag(s) <- 1
  s <- s * spread
  mu <- mvtnorm::dmvnorm(coordinates[, 1:3], mean = location,
                         sigma = s)
  mu <- mu * targetSnr/max(mu)
  coordinates$signal <- round(mu, 3)

  # Generating noise + data -----------------
  pop_n = dim(dat)[4]
  samp = sample(pop_n, grp_size*2, replace = FALSE)
  noise <- as.numeric(apply(dat[,,,samp[1:grp_size]],c(1,2,3),mean) -
                         apply(dat[,,,samp[grp_size+(1:grp_size)]],c(1,2,3),mean))

  scale_noise = sd(noise)
  coordinates$noise <- noise/scale_noise
  coordinates$scale_coef = scale_noise
  coordinates$observed <- coordinates$signal + coordinates$noise

  return(coordinates)
}
