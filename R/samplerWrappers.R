
# This file has service functions that wrap different truncated mutlivariate gaussian
# samplers
# Sampling functions can be added, provided they use the same arguments. 
# Arguments not needed can be handled using '...'
#




setSampler = function(sampler){
  if (sampler == "restrictedMVN") {     #Recomended sampler!!
    samp_func = wrap_restrictedMVN
  } else if (sampler == "rtmvnorm") {
    samp_func = wrap_rtmvnorm
  } else {
    stop("Please specify sampler ('restrictedMVN' / 'hamiltonian' / 'rtmvnorm')")
  }
  return(samp_func)
}

wrap_restrictedMVN = function(n, mu, sigma, lower = rep( -Inf, length = length(mu)),
                          upper = rep( Inf, length = length(mu)), start.value = NULL,
                          burn.in = length(mu)*1000, thinning = 1, algorithm = NULL,  # just to catch this argument
                          ...) {

  constr = restrictedMVN::thresh2constraints(d = length(mu),
                                              lower = lower, upper = upper)

  samp = restrictedMVN::sample_from_constraints(linear_part = constr$linear_part,
                                 offset= constr$offset,
                                 mean_param = mu,
                                 covariance = sigma,
                                 initial_point = start.value,                                 ndraw=n*thinning,
                                 burnin=burn.in)
  thin_samp = samp[(1:n)*thinning,]
  return(thin_samp)
}


####
#
# Wrapper for a simple rejection sampler 'tmvtnorm' that is quite stable for small regions
#
#
####

wrap_rtmvnorm = function(n, mu, sigma,lower = rep( -Inf, length = length(mu)),
                         upper = rep( Inf, length = length(mu)), ...) {
  # Sampling based on the quite stable tmvtnorm by Wilhelm and Manjunath
  # For stable sampling of small examples, this is recomended.
  return(tmvtnrom::rtmvnorm(n, sigma = sigma, mean = mu, lower = lower,upper = upper, ...))
}

