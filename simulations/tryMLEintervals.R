
require('parallel')

thresh = 0.08   # threshold
grp_size = c(4,8,16)  # group size (and variance)
mu_size = 0:20 *0.02 # Set avg difference between groups

# Set shape of external differences
mu_shape = matrix(nr = 2, c(c(0.5,rep(1,3), 0.5), c(0,rep(1,3), 0)))

# Set covariance
sigCov = list()
sigCov[[1]] = matrix(c(0.04, 0.02, 0.006, 0.0, 0.0,
                       0.02,  0.04, 0.016, 0.0, 0.0,
                       0.006,  0.016,  0.03, 0, 0,
                       0., 0.0,  0, 0.04, 0.01,
                       0.0, 0.0, 0 , 0.01 , 0.03),nrow= 5)
sigCov[[2]] = diag(5)*0.04

# Set whether covariance is known
sigKnown = c(TRUE)


ntests = 10
n_musize = length(mu_size)
n_grps = length(grp_size)
n_shapes = nrow(mu_shape)
n_sigs = length(sigCov)
n_known = length(sigKnown)

n_exp = prod(c(n_musize,n_grps,n_shapes,n_sigs,n_known))

samples = list()
confIntsExact = list()
confIntsMLE = list()

i_exp = 0

for (grp_i in 1:n_grps){
  for (musize_i in 1:n_musize){
    for (shape_i in 1:n_shapes){
      for (sig_i in 1:n_sigs){
        # Parameters needed to generate the sample...
        dat_ind = (grp_i-1)*n_musize*n_shapes*n_sigs+
          ((musize_i-1)*n_shapes*n_sigs)+ (shape_i-1)*n_shapes+sig_i
        mu_vec = mu_shape[shape_i,]*mu_size[musize_i]
        sampParams = setSamplingParameters(na=grp_size[grp_i], nb = grp_size[grp_i],
                                           mudiff=mu_vec,
                                           sig = sigCov[[sig_i]],
                                           conditioning = c(-1,1,1,1,-1))
        # Generate Data
        set.seed(dat_ind)
        samples[[dat_ind]] = mclapply(rep(1,ntests),condSample,
                                      sP = sampParams, thresh = thresh,
                                      returndat = TRUE, mc.cores = useCores)

        # Make inference
        for (known_i in n_known) {
          useSig = NULL
          if (sigKnown[known_i]==1) {
            useSig = diffCov(sampParams$sig, samples[[dat_ind]][[1]]$groups)
          }
          i_exp = i_exp + 1

          confIntsMLE[[i_exp]] = list()
          for(i in 1:ntests){
            confIntsMLE[[i_exp]][[i]] <-
              optimizeSelected(as.numeric(samples[[dat_ind]][[i]]$condZ),
                                     useSig, thresh,
                                     selected = c(FALSE,TRUE,TRUE,TRUE,FALSE),
                                     stepRate = 0.6,
                                     delay = 100,
                                     assumeConvergence = 3000,
                                     trimSample = 100,
                                     maxiter = 7000)
          }
          ### Here is where we put the inference function
          confInts[[i_exp]] = mclapply(1:ntests, confFunc,
                                                   samp = samples_b[[dat_ind]],
                                                   sampParams = sampParams,
                                                   Sigma = useSig,
                                                   cond_vec = condExternal[1,],
                                                   regularize = regularize,
                                                   mc.cores = useCores)

        }
        cat(".")
      }
    }
  }
}
