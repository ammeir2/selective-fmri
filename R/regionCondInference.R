
# Main source file for inference on regions

#source('samplerWrappers.R')


regionConfInt = function(dat, groups, threshold, Sigma = NULL, 
                      CIparams, conditioning, external_mu_val=NULL,
                      tilt = FALSE,
                      regularize_sig = 0,
                      sampler = "restrictedMVN",...){
  # We assume the parameter exceeds the threshold, 
  # and find a confidence interval for the sum
  
  datSt = getDiffMoments(as.matrix(dat), groups , regularize_sig)
  
  if (! is.null(Sigma)) {
    datSt$cov_diff = Sigma
    CIparams$sig_known = TRUE
  } else{ 
    CIparams$sig_known = FALSE
    
    
  }
  
  CI_res = formTwosidedCI(diffVec = datSt$delta_hat,thresh = threshold,  conditioning = conditioning,                                                  
                            diffSigma = datSt$cov_diff,external_mu_val = external_mu_val,
                            tilt = tilt,params =CIparams, sampler = sampler,...)                


  return(CI_res)
}


formTwosidedCI = function(diffVec,thresh,diffSigma,  params,
                            sampler = "hamiltonian",conditioning = rep(1,length(diffVec)),
                            inner_shape = NULL,
                            external_mu_val=0,
                            start_search = "obs",
                            tilt = FALSE,
                            median_sweep = FALSE,
                            doPlot = FALSE){
  
  #
  # thresh is the thresholding value (always positive). 
  # A conditioning vector of {-1,0,1}^d is given in the input. 
  
  # For a positive region:
  #   condition[i] = 1 means diffVec[i]>thresh
  #   condition[i] = -1 means diffVec[i]<thresh
  
  # For a negative region:
  #   condition[i] = 1 means diffVec[i] < -thresh
  #   condition[i] = -1 means diffVec[i] > -thresh
  
  # eta is taken to be 1 if condition[i]=1 and 0 otherwise
  # An unchecked assumption is that the region 'conditioning==1' is a single connected region.
  
  thin_c = 6 # thinning per variable
  
  ###
  # 1. Preparation
  ###
  alpha = params$alpha
  d = length(diffVec)
  eta = 1*(conditioning==1)
  region = which(conditioning==1)
  exterior = which(conditioning==-1)
  
  diffSigmaInv = solve(diffSigma)
  
  stopifnot(sum(eta)>0)  # at least one positive point in the region
  
  # Choose sampler function
  samp_func = setSampler(sampler)
  
  # The inference function runs for positive regions, so we first set direction
  # Check direction
  direction = 1
  if (diffVec[region[1]] < 0) {
    diffVec = -diffVec
    direction = -1
  }
  stopifnot(all ( diffVec[region]>= thresh) & all(diffVec[exterior] < thresh))
  
  # Set lower and upper thresholds for conditioning
  lower = rep(-Inf, d)
  lower[region] = thresh
  
  upper = rep(Inf, d)
  upper[exterior] = thresh
  
  # Find the observed statistic (area-of-region)
  stat_obs = sum(eta * diffVec)
  
  obs_quant = numeric(length(params$m_values))
  
  # make sure the profile has a sum of 1. 
  if (is.null(inner_shape)) {
    inner_shape = getMuSigmaEta(1, diffSigma[region,region], eta[region])
  }
  inner_shape = inner_shape/mean(inner_shape)

  if (is.null(external_mu_val)){
    external_mu_val = diffVec[exterior]
  } 
  

  # Choose starting value:  
  if (is.numeric(start_search)){
    m0 = start_search
  } else if (start_search == 'thresh'){
    m0 = (thresh + stat_obs/sum(eta))/2
  } else if (start_search == 'obs') {
    m0 = stat_obs/sum(eta)
  }

  start.value = diffVec
  
  if (tilt) {
    
    # We want to search for a good starting value of m0, so stat_obs is the sum of the region
    
    m0stop = FALSE
    m0counter = 0  
    
    full_mus = numeric(length(eta))
    while(m0stop){  
      
      mean_vec = setMeanVec(m0, inner_shape, external_mu_val, conditioning)
      search_samp = samp_func(50, mu = mean_vec, sigma = diffSigma, lower = lower, upper = upper,
                              burn.in = d*100*thin_c,thinning = d*thin_c,algorithm = "gibbs",
                              start.value = start.value)
    
      # Check that the observed statistic is not in the tails of the distribution for this m0
      search_res = searchStep( search_samp, eta, stat_obs, quants = c(0.2,0.8), m0, params$searchstep) 

      # Check if failed to get a non NA sample
      if(is.na(search_res$success)){  return(list())  }
      
      # process search results; stop if the following conditions occur:
      # (m0 < min(params$m_values))  | (m0 > max(params$m_values) | (m0counter > 20)  )
      m0counter = m0counter + 1
      if ((m0 < min(params$m_values))  | (m0 > max(params$m_values) | (m0counter > 20)  )){
        m0stop = TRUE
        force_break = TRUE
      } else{
        m0stop = search_res$success
        m0 = search_res$new_m0
      }
    }

    mean_vec = setMeanVec(m0, inner_shape, external_mu_val, conditioning)
    samp = samp_func(params$nsamp, mu = mean_vec, sigma = diffSigma, lower = lower, upper = upper,  
                     burn.in = d*100*thin_c, thinning = d*thin_c , algorithm = "gibbs", 
                     start.value = start.value)
    # Check if failed to get a non NA sample
    if(any(is.na(samp))){  return(list())  }


    ts = samp %*% eta
    o_samp = samp[order(ts),]
    o_ts = o_samp %*% eta
    
    mu_shape = numeric(length(eta))
    mu_shape[region] = inner_shape
    
    weights = o_samp %*% t(mu_shape %*% diffSigmaInv) 
    # Note that when weights are too small for tilting, sometimes we might get Nan's
  
    for (i in 1:length(params$m_values)){
      curr_m = params$m_values[i]
      min_wt = min((curr_m-m0)*weights)
      rescaled_weights = exp((curr_m-m0)*weights-min_wt)
      norm_weights = rescaled_weights/sum(rescaled_weights)
      obs_quant[i] = sum(norm_weights[o_ts<stat_obs])
    }
    
  } else {
    full_mus = numeric(length(eta))
    for (i in 1:length(params$m_values)){
      mean_vec = setMeanVec(params$m_values[i], inner_shape, external_mu_val, conditioning)

      samp = samp_func(params$nsamp, mu = mean_vec, sigma = diffSigma, lower = lower, upper = upper,
                         burn.in = d*100*thin_c,thinning = d*thin_c,algorithm = "gibbs",
                         start.value = diffVec)
      if (any(is.na(samp))){      
        samp = samp_func(params$nsamp, mu = mean_vec, sigma = diffSigma, lower = lower, upper = upper,
                           burn.in = d*100*thin_c,thinning = d*thin_c,algorithm = "gibbs",
                           start.value = (diffVec + thresh)/2)
      }    
    
      ts = samp %*% eta
      # find quantile of stat_obs in o_ts:
      obs_quant[i] = (sum(stat_obs >= ts)+1)/(length(ts)+2)
    }
  }  

  # sanity check
  if(!any(obs_quant<0.98 & obs_quant>0.02,na.rm = TRUE)){
    return(list())
    # Failed run
  }
  
  index_0 = which(params$m_values==0)
  p_value = 1 - obs_quant[index_0]
  
  lower_bound = findMValue(params$m_values,obs_quant,1-alpha/2)
  upper_bound = findMValue(params$m_values,obs_quant,alpha/2)
  reject_reg = (params$m_values < lower_bound) | (params$m_values > upper_bound)
  
  mean_est= findMValue(params$m_values,obs_quant,0.5)

  full_baseline = setMeanVec(1, inner_shape,  external_mu_val, conditioning)
  
  return( list(vals = params$m_values, rej = reject_reg, quants = obs_quant, 
               baseline = full_baseline, params = params, 
               direction = direction,p_value = p_value, stat_obs = stat_obs,
               lower_bound = lower_bound, upper_bound = upper_bound, cond_mean_est = mean_est))
}



getDiffMoments = function(datmat, groups, regularize_sig=0) {
# Extracts the for a difference-of-means model information from the samples, 
# including :
# - Vector of estimators Z
# - Estimated covariance per group
# - Estimated covariance of difference

  reg_cov_weight = regularize_sig
  
  datSt = list()
  if(nrow(datmat)==1){
    datmat = t(datmat)
  } 
  stopifnot(nrow(datmat)==length(groups))
  
  datSt$n_a = sum(groups==1)
  datSt$n_b = sum(groups==0)

  if (ncol(datmat)>1) {
    dat_a = datmat[groups==1,]
    dat_b = datmat[groups==0,]
    
    datSt$cov_a = cov(dat_a)
    datSt$cov_b = cov(dat_b)
    
  # to regularize sigma, we inflate by (1+reg_cov_weight) 
    reg_cov_a = datSt$cov_a
    diag(reg_cov_a) = diag(reg_cov_a)*(1+reg_cov_weight)
    datSt$cov_a = reg_cov_a
      
    reg_cov_b  = datSt$cov_b
    diag(reg_cov_b) = diag(reg_cov_b)*(1+reg_cov_weight)
    datSt$cov_b = reg_cov_b
      
  }
  else{
    dat_a = datmat[groups==1]
    dat_b = datmat[groups==0]
    
    datSt$cov_a = var(dat_a)*(1+reg_cov_weight)
    datSt$cov_b = var(dat_b)*(1+reg_cov_weight)
    
    datSt$delta_hat = mean(dat_a) - mean(dat_b)
  }
  datSt$cov_diff = datSt$cov_a/datSt$n_a + datSt$cov_b/datSt$n_b
  datSt$delta_hat = colMeans(dat_a) - colMeans(dat_b)
  
  return(datSt)    
}

setMeanVec = function(m, inner_shape, ext_shape, conditioning){
  # Takes the profile vector inner_shape and scales by m
  # Then adds the ext_shape to the external region
  mean_vec = numeric(length(conditioning))
  region = which(conditioning==1)
  mean_vec[region] = m *inner_shape
  
  exterior = which(conditioning==-1)
  mean_vec[exterior] = ext_shape 
  return(mean_vec)
}

searchStep = function ( samp, eta, stat_obs, quants = c(0.2,0.8), m0, step_size) {
  # Takes a sample and decides whether m0 is in the qunatile range defined 
  # by quants, or needs to be changed (by step_size)
  
  stopifnot(length(quants)==2)
  
  # check for NAs
  search_res = list()
  
  if(any(is.na(samp))){
    # Failed to get a non NA sample
    search_res$success = NA
    search_res$new_m0 = m0
    return(search_res)
  }
  
  search_ts = samp %*% eta
  search_qs = quantile(ecdf(search_ts),c(quants))
  
  if ((stat_obs < search_qs[2] ) & (stat_obs > search_qs[1] )) {
    search_res$success = TRUE
    search_res$new_m0 = m0
  }
  else {
    search_res$success = FALSE
    if (stat_obs > search_qs[2]) { 
      search_res$new_m0 = m0 +step_size
    }  else {
      search_res$new_m0 = m0 - step_size
    } 
  }
  return(search_res)
}



getMuSigmaEta = function(m, Sigma,eta) {
  # Get the mean vector from the Sigma and a parameter m0. 
  # The goal is to get a positive mean vector mu = m*shape
  # so that m is the magnitude:
  #     m * sum(support(eta)) = eta'mu    ->   m  =  eta'mu/support(eta)
  # and its shape is proportion to eta'Sigma 
  #    shape propto eta'Sigma
  
  if (is.matrix(Sigma)){
    shape =  as.vector(eta %*% Sigma)
    normalizing = sum(eta * shape)/sum(eta)
    shape = shape/normalizing
  } else {  # in case Sigma is a real-value
    shape = 1
  }
  return(m * shape)
}


findMValue = function(m_values,obs_quant,p){
  # A bit counter intuitive, but as m increases the reference distribution shifts right.
  # Therefore, the lower bound for m is should be the first value where obs_quant is smaller than p,
  # and the upper bound for m is the first value where obs_quant is larger than p. 
  
  # If we use no tilt, we might not have monotone values in the quantile estimates. 
  # Instead of smoothing, we give a score to how many values agree with this selection. 
  
  n = length(obs_quant)
  score = numeric(n+1)
  score[1] = sum(obs_quant< p, na.rm=TRUE)
  for (i in 2:length(score)){
    score[i] = sum(obs_quant[1:(i-1)] >= p,na.rm=TRUE) + sum(obs_quant[(i):n]< p,
                                                                      na.rm=TRUE)
  }
  score[n+1] = sum(obs_quant > p,na.rm=TRUE)
  best_match = which.max(score)
  quantile = c(-Inf,m_values)[best_match]
  return(quantile)
}

