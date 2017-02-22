
#### 
#
#sampParams = setSamplingParameters(na=19, nb = 17, sig = diag(3) *0.04, mudiff=0.3)
#
####

setSamplingParameters = function(na,sig, nb = NULL, sigb = NULL, mudiff, conditioning = 1) {
  # mudiff is mean difference between group a and group b
  # we set the mean of group b to 0, and group a is mudiff

  # Fill in group b with defaults 
  if (is.null(nb)) { nb = na}
  if (is.null(sigb)) { sigb = sig }
  stopifnot(ncol(sig)==ncol(sigb))

  # expand mudiff
  d = ncol(sig)
  if (length(mudiff) == 1){ mudiff = rep(mudiff,d) }
  stopifnot(length(mudiff) == d )
  
  # set param structure
  sampParam = list(na = na, nb = nb, sig = sig, sigb = sigb, mudiff = mudiff,d = ncol(sig), 
                   conditioning = conditioning)
  return(sampParam)
}

#B = 500   # Parametric bootstrap samples

# sampP = setSamplingParameters(na=20, sig = diag(3) *0.04, mudiff=0.3)
# condSample(B = 500, sP = sampP,thresh=0.5)
condSample = function(B,sP, thresh, returndat = FALSE){
  
  # The Y matrix is the full data matrix
  # Y_a are samples for group a, Y_b for group b
  # Z is the mean-difference vector
  # cond means conditional on thresholding. 
  # sP$conditioning is a vector: conditioning[i] = 1 means Z[i]>thresh
  #                           conditioning[i] = -1 means Z[i] < thresh
  #                           conditioning[i] = 0 means Z[i] no condition
  #              it is checked using all((conditioning * Z) >= (conditioning * thresh)
  
  # Simulate the conditional distribution 

  groups = rep(c(1,0),times = c(sP$na,sP$nb))

  condZ = matrix(nr = B, nc = sP$d)
  if (returndat){
    condY = array(0, c(sP$na+sP$nb,sP$d,B))
  } else{
    condY = NULL
  }
  groups = c(rep(1,sP$na),rep(0,sP$nb))

  if (is.null( sP$conditioning)){
    conditioning = 1
  } else {
    conditioning = sP$conditioning
  }
  stopifnot(is.element(length(conditioning), c(1, sP$d)))
  
  B_i = 0   # counter passed conditioning
  i = 0     # global counter
  while(B_i<B){
    i = i + 1
    Y_a = mvrnorm(n = sP$na,mu = sP$mudiff,Sigma = sP$sig)
    Y_b = mvrnorm(n = sP$nb,mu = rep(0, sP$d),Sigma = sP$sigb)
    Z = colMeans(Y_a) - colMeans(Y_b)

    # check condition:
    check_cond = all(Z*conditioning >= thresh * conditioning)
    if (check_cond) {
      B_i = B_i + 1
      condZ[B_i,] = Z 
      if (returndat){
        condY[,,B_i] = rbind(Y_a, Y_b)
      }
    }
  }
  return(list(condZ = condZ, samples =  i, condY= condY, groups = groups))
}


diffCov = function(sampleCov, groups) {
  n_a = sum(groups==0)
  n_b = sum(groups==1)
  return(sampleCov * (1/n_a + 1/n_b) )
}



