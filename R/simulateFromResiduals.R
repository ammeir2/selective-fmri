# #library(selectivefmri)
#
# residualData3D <- function(dat, group_ns, targetSnr ,spread = 1) {
#
#     dims = dim(dat)
#
#     # Generating signal
#     coordinates <- expand.grid(i=1:dims[1],j=1:dims[2],k=1:dims[3])
#     nnodes <- prod(dims[1:3])
#     coordinates$row = 1:nrow(coordinates)
#
#     location <- sapply(dims[1:3], function(x) sample.int(x,1))
#     s <- matrix(0.3, nrow = 3, ncol = 3)
#     diag(s) <- 1
#     s <- s * spread
#     mu <- mvtnorm::dmvnorm(coordinates[, 1:3], mean = location,
#                            sigma = s)
#     mu <- mu * targetSnr/max(mu)
#     coordinates$signal <- round(mu, 3)
#
#     # Generating noise + data -----------------
#     pop_n = dim(dat)[4]
#     samp = sample(pop_n, group_ns*2, replace = FALSE)
#     noise <- as.numeric((apply(dat[,,,samp[1:group_ns]],c(1,2,3),mean) -
#         apply(dat[,,,samp[group_ns+(1:group_ns)]],c(1,2,3),mean))/2)
#
#     noise_try = matrix(nc =2,nr=1000)
#     for (i in 1:1000){
#       ss = sample(pop_n, group_ns*2, replace = FALSE)
#       noise_try[i,] = as.numeric((apply(dat[,,,ss[1:group_ns]],c(1,2,3),mean) -
#                     apply(dat[,,,ss[group_ns+(1:group_ns)]],c(1,2,3),mean))/2)[1:2]
#     }
#     scale_noise = sd(noise)
#     coordinates$noise <- noise/scale_noise
#     coordinates$scale_coef = scale_noise
#     coordinates$observed <- coordinates$signal + coordinates$noise
#
#     return(coordinates)
# }
#
#
# run.sim.cube <- function(config, dat) {
#   snr <- as.numeric(config["snr"])
#   spread <- as.numeric(config["spread"])
#   p_threshold <- as.numeric(config["BHlevel"])
#   replications <- as.numeric(config["replications"])
#   group_ns <-as.numeric(config["ns"])
#
#   dims <- dim(dat)[1:3]
#   pop_size <- dim(dat)[4]
#
#   # Estimate the covariance
#
#   dim(dat) <- c(prod(dims),pop_size)
#   tdat = t(dat)
#   alpha <- 1/pop_size
#   sds = apply(tdat, 2, sd)
#   sample_covariance <- (1-alpha)*cov(tdat)+alpha*diag(sds^2)
#   covariance <- sample_covariance*(1/group_ns*2)
#   dim(dat) <- c(dims[1:3],pop_size)
#   rm(sample_covariance)
#   rm(tdat)
#
#   # Start the main loop
#   simresults <- list()
#   simcover <- matrix(nrow = replications, ncol = 2)
#   colnames(simcover) <- c("cover", "power")
#   for(rep in 1:replications) {
#     selected <- FALSE
#     maxsize <- 0
#     while(is.na(maxsize) | maxsize < 2) {
#       coordinates <- residualData3D(dat, group_ns, snr, spread)
#       coordinates$zval <- (coordinates$observed / sds) *
#                                   sqrt(group_ns) * coordinates$scale_coef
#       coordinates$pval <- 2 * pnorm(-abs(coordinates$zval))
# #      coordinates$qval <- p.adjust(coordinates$pval, method = "BH")
#       coordinates$selected <- coordinates$pval < p_threshold
#       selected <- coordinates$selected
#       if(sum(selected) >= 2) {
#         clusters <- findClusters(coordinates)
#         maxsize <- max(sapply(clusters, function(x) sum(x$selected)))
#       }
#       print('.')
#     }
#
#     print(c(rep = rep,  p_threshold = p_threshold, snr = snr, spread = spread))
#     print(c(nselected = sum(selected)))
#     sizes <- sapply(clusters, nrow)
#     threshold <- abs(qnorm(p_threshold))
#
#     results <- list()
#     iterCover <- 0
#     naiveCover <- 0
#     profCover <- 0
#     weights <- 0
#     slot <- 1
#
#     mse <- 0
#     msecount <- 0
#     for(m in 1:length(clusters)) {
#       #results[[m]] <- list()      ## Erase line?
#       cluster <- clusters[[m]]
#       cluster <- subset(cluster, !is.na(cluster$selected))
#       selected <- coordinates$selected[cluster$row]
#       observed <- coordinates$observed[cluster$row]
#       signal <- coordinates$signal[cluster$row]
#       subCov <- covariance[cluster$row, cluster$row, drop = FALSE]
#
#
#       print(c(round(m / length(clusters), 2), nrow(cluster)))
#
#       try(mle <- optimizeSelected(observed, subCov, threshold,
#                                   selected = selected,
#                                   projected = NULL,
#                                   stepRate = 0.6,
#                                   coordinates = cluster[, 1:3],
#                                   tykohonovParam = NULL,
#                                   tykohonovSlack = 2,
#                                   stepSizeCoef = 2,
#                                   delay = 10,
#                                   assumeConvergence = 1000,
#                                   trimSample = 200,
#                                   maxiter = 1200,
#                                   probMethod = "selected",
#                                   init = observed,
#                                   imputeBoundary = "neighbors"))
#
#       conditional <- mean(mle$conditional[selected])
#       naive <- mean(observed[selected])
#       true <- mean(signal[selected])
#       mse <- mse * msecount / (msecount + 1) + c(naive = (naive - true)^2, conditional = (conditional - true)^2) / (msecount + 1)
#       msecount <- msecount + 1
#
#
#       for(methodind in 1:2) {
#         if(methodind == 1) {
#           method <- "onesided"
#         } else if(methodind == 2) {
#           method <- "selected"
#         }
#         results[[slot]] <- list()
#         print(c(round(m / length(clusters), 2), nrow(cluster)))
#         subCov <- covariance[cluster$row, cluster$row, drop = FALSE]
#
#         # Computing the p-value based on samples from the null
#         nullfit <- NULL
#         try(nullfit <- optimizeSelected(observed, subCov, threshold,
#                                         projected = 0,
#                                         selected = selected,
#                                         stepRate = 0.6,
#                                         coordinates = cluster[, 1:3],
#                                         tykohonovParam = NULL,
#                                         tykohonovSlack = 0.0001,
#                                         stepSizeCoef = 0,
#                                         delay = 1,
#                                         assumeConvergence = 1,
#                                         trimSample = 200,
#                                         maxiter = 1500,
#                                         probMethod = method,
#                                         init = rep(0, length(observed)),
#                                         imputeBoundary = "neighbors"))
#
#
#         if(is.null(nullfit)) next
#         naive <- mean(observed[selected])
#         samp <- nullfit$sample
#         nullmeans <- NULL
#         try(nullmeans <- rowMeans(samp[, selected, drop = FALSE]))
#         if(is.null(nullmeans)) next
#         # p-value based on one sided test
#         if(naive > 0) {
#           pvalue <- mean(naive < nullmeans)
#         } else {
#           pvalue <- mean(naive > nullmeans)
#         }
#
# #        obsratio <- abs(mean(observed[selected]) / mean(signal[selected]))
#
#         try(profile <- optimizeSelected(observed, subCov, threshold,
#                                         selected = selected,
#                                         projected = mean(signal[selected]),
#                                         stepRate = 0.6,
#                                         coordinates = cluster[, 1:3],
#                                         tykohonovParam = NULL,
#                                         tykohonovSlack = 2,
#                                         stepSizeCoef = 2,
#                                         delay = 10,
#                                         assumeConvergence = 750,
#                                         trimSample = 200,
#                                         maxiter = 1800,
#                                         probMethod = method,
#                                         init = observed,
#                                         imputeBoundary = "neighbors"))
#
#         samp <- profile$sample
#         profMeans <- NULL
#         try(profMeans <- rowMeans(samp[, selected, drop = FALSE]))
#         if(is.null(profMeans)) next
#         profPval <- 2* min(mean(naive < profMeans), mean(naive > profMeans))
#
#         true <- mean(signal[selected])
#         profResult <- c(true = mean(signal[selected]), profPval = profPval, pvalue = pvalue)
#         results[[slot]][[1]] <- c(snr = snr, spread = spread, method = methodind,p_threshold = p_threshold,
#                                   group_ns = group_ns,size = sum(selected), true = true)
#         results[[slot]][[2]] <- profResult
#         results[[slot]][[3]] <- c(conditional = conditional, naive = naive, true = true)
#
#         print(profResult)
#
#         weight <- 1
#         iterCover <- iterCover + weight * (pvalue < 0.05)
#         profCover <- profCover + weight * (profPval > 0.05)
#         weights <- weights + weight
#
#         slot <- slot + 1
#
#       }
#     }
#
#     simcover[rep, 1] <- sum(profCover) / sum(weights)
#     simcover[rep, 2] <- sum(iterCover) / sum(weights)
#     print(c(iterprof = profCover, itermle = iterCover) / weights)
#     print(colMeans(simcover[1:rep, , drop = FALSE]))
#     simresults[[rep]] <- results
#   }
#   return(simresults)
# }
#
# #load('/Users/yuvalb/Dropbox/SelectiveFmri/brain_data_4mm_Cambridge.rda')
# #configurations <- expand.grid(snr = c(0.8, 0.05, 0.2),BHlevel = 0.1,replications = 50,ns = c(8,16,32),spread=2)
#
# #set.seed(500)
# #system.time(simResults <- apply(configurations, 1, run.sim.cube,brain_data$coef_cube))
# #save(simResults, file = "simulations/results/Resid_Apr.Robj")
#
#
# process_results = function(res){
#   # Processing ------------------------
# #  require(ggplot2)
# #  require(dplyr)
# #  require(reshape2)
#   resultList <- list()
#   len <- 1
#   for(i in 1:length(res)) {
#     for(j in 1:length(res[[i]])) {
#       nonEmpty <- sapply(res[[i]][[j]], length) > 0
#       iterResult <- data.frame(matrix(nrow = sum(nonEmpty), ncol = 13))
#       iterList <- res[[i]][[j]]
#       row <- 1
#       exp <- j * runif(1)
#       for(k in which(nonEmpty)) {
#         method <- iterList[[k]][[1]][[3]]
#         if(method == 1) {
#           method <- "onesided"
#         } else if(method == 2) {
#           method <- "selected"
#         }
#         iterResult[row, ] <- c(experiment = exp,
#                                cluster = row,
#                                snr = as.numeric(iterList[[k]][[1]][["snr"]]),
#                                spread = as.numeric(iterList[[k]][[1]][["spread"]]),
#                                method = method,
#                                rho = -1,
#                                BHlevel = as.numeric(iterList[[k]][[1]][["BHlevel"]]),
#                                ns = as.numeric(iterList[[k]][[1]][["ns"]]),
#                                size = as.numeric(iterList[[k]][[1]][["size"]]),
#                                true = as.numeric(iterList[[k]][[1]][["true"]]),
#                                pval = as.numeric(iterList[[k]][[2]][["pvalue"]]),
#                                cover = as.numeric(iterList[[k]][[2]][["profPval"]]),
#                                naive = iterList[[k]][[3]][["naive"]],
#                                cond = iterList[[k]][[3]][["conditional"]])
#         row <- row + 1
#       }
#       iterResult <- data.frame(iterResult)
#       names(iterResult) <- c("experiment", "cluster", "snr", "spread", "method", "rho",
#                              "pthreshold", "size", "true", "pval", "cover", "naive", "cond")
#       resultList[[len]] <- iterResult
#       len <- len + 1
#     }
#   }
# }
