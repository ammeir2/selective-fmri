set.seed(505)
snr <- 4.5
spread <- 2
BHlevel <- 0.01
slack <- 1
noise_type <- "fmri"
config <- list()
config[["grp_size"]] <- 16

load('fmridata/brain_data_4mm_Cambridge.rda')
noise_dat <- brain_data$t_cube[1:11,1:11,1:9,]


if(noise_type == "sim" ) {
  rho <- config[["rho"]]
  grp_size <- -1
} else {
  grp_size <- config[["grp_size"]]
  rho <- -1
}

I <- 11
J <- 11
K <- 9
dims <- c(I, J, K)

coordinates <- expand.grid(i = 1:I, j = 1:J, k = 1:K)

if (noise_type == "sim") {
  covariance <- rho^as.matrix(dist(coordinates[, 1:3], method = "euclidean",
                                   diag = TRUE, upper = TRUE))
  covEigen <- eigen(covariance)
  sqrtCov <- covEigen$vectors %*% diag(sqrt(covEigen$values)) %*% t(covEigen$vectors)
} else if (noise_type == "fmri") {
  stopifnot(all.equal(dim(noise_dat)[1:3], dims))
  pop_size <- dim(noise_dat)[4]
  dim(noise_dat) = c(prod(dims), pop_size)
  alpha <- 1/pop_size
  sample_sds <- apply(noise_dat, 1, sd)
  # covariance of a single sample
  sample_covariance <- (1-alpha)*cov(t(noise_dat))+alpha*diag(sample_sds^2)
  # covariance of group averages
  covariance <- sample_covariance*(1/grp_size*2)
  dim(noise_dat) <- c(dims[1:3],pop_size)
  rm(sample_covariance)
}

simresults <- list()
origCov <- covariance
covariance <- origCov
maxsize <- 0
while(is.na(maxsize) | maxsize < 2) {
  if (noise_type == "sim") {
    coordinates <- generateArrayData3D(dims, sqrtCov, snr, spread)
    sds <- 1
  }
  else if (noise_type == "fmri"){
    coordinates <- residualData3D(noise_dat, grp_size, snr, spread)
    sds <- sqrt(diag(covariance))/coordinates$scale_coef
  }

  coordinates$observed <- coordinates$observed / sds
  coordinates$signal <- coordinates$signal / sds
  coordinates$zval <- coordinates$observed
  coordinates$pval <- 2 * pnorm(-abs(coordinates$zval))
  coordinates$qval <- p.adjust(coordinates$pval, method = "bonferroni")
  coordinates$qval <- coordinates$pval  # FOR NO CORRECTION!!!!
  coordinates$selected <- coordinates$qval < BHlevel
  selected <- coordinates$selected
  if(sum(selected) >= 2) {
    threshold <- abs(qnorm(BHlevel / 2))
    coordinates$selected <- coordinates$zval > threshold
    pclusters <- findClusters(coordinates)

    coordinates$selected <- coordinates$zval < -threshold
    nclusters <- findClusters(coordinates)
    clusters <- c(pclusters, nclusters)
    coordinates$selected <- coordinates$qval < BHlevel

    maxsize <- max(sapply(clusters, function(x) sum(x$selected)))
  }
}

covariance <- cov2cor(covariance)

print(c(rep = rep, rho = rho, grp_size = grp_size, BHlevel = BHlevel, snr = snr, spread = spread))
print(c(nselected = sum(selected)))
sizes <- sapply(clusters, nrow)

results <- list()
iterPower <- 0
iterCover <- 0
weights <- 0
slot <- 1
mse <- 0
msecount <- 0
resultList <- list()
totalSelected <- sum(sapply(clusters, function(x) sum(x$selected)))
for(m in 1:length(clusters)) {
  cluster <- clusters[[m]]
  cluster <- subset(cluster, !is.na(cluster$selected))
  selected <- coordinates$selected[cluster$row]
  observed <- coordinates$observed[cluster$row]
  signal <- coordinates$signal[cluster$row]
  subCov <- covariance[cluster$row, cluster$row, drop = FALSE]
  naive <- mean(observed[selected])

  print(c(round(m / length(clusters), 2), nrow(cluster)))
  alpha <- 0.1 * sum(selected) / totalSelected

  fit <- selectiveMRI(observed, subCov, threshold, contrast = NULL,
                      coordinates = cluster[, 1:3],
                      selected = selected,
                      tykohonovSlack = 1,
                      CIalpha = alpha,
                      probMethod = "selected",
                      imputeBoundary =  "neighbors",
                      control = selectiveMRI_control())

  e <- rep(0, length(selected))
  e[selected] <- rep(1 / sum(selected), sum(selected))
  naivesd <- sqrt(as.numeric(t(e) %*% subCov %*% e))
  naiveCI <- naive + c(-1, 1) * qnorm(1 - alpha / 2)

  fit$true <- mean(signal[selected])
  fit$size <- sum(selected)
  fit$naive <- naive
  fit$conditionSize <- length(selected) - sum(selected)
  fit$naiveCI <- naiveCI
  resultList[[m]] <- fit
}

#saveRDS(fit, "misc scripts/examples/ci example seed 505.rds")

conditional <- sapply(resultList, function(x) x$conditional)
true <- sapply(resultList, function(x) x$true)
naive <- sapply(resultList, function(x) x$naive)
cbind(conditional, naive, true)
CI <- t(sapply(resultList, function(x) x$CI))
cbind(CI[, 1], true, CI[, 2])
naiveCI <- t(sapply(resultList, function(x) x$naiveCI))
size <- sapply(resultList, function(x) x$size)

# Plotting --------------------
offset <- 0.2
naivedat <- data.frame(estimate = naive, nselected = size,
                       var  = 1:length(naive),
                       lCI = naiveCI[, 1], uCI = naiveCI[, 2],
                       method = "naive", offset = offset)
conddat <- data.frame(estimate = conditional, nselected = size,
                      var  = 1:length(naive),
                      lCI = CI[, 1], uCI = CI[, 2],
                      method = "mle", offset = 0)
truedat <- data.frame(estimate = true, nselected = size,
                      var  = 1:length(naive),
                      lCI = NA, uCI = NA,
                      method = "true", offset = offset * 2)
forplot <- rbind(naivedat, conddat, truedat)
forplot$method <- factor(forplot$method, levels = c("mle", "naive", "true"))

ggplot(forplot) +
  geom_point(aes(x = var + offset, y = estimate, shape = method, size = nselected)) +
  geom_segment(data = subset(forplot, method != "true"), aes(x = var + offset, xend = var + offset,
                   y = lCI, yend = uCI, col = method, linetype = method)) +
  theme_bw() + geom_hline(yintercept = 0) + xlab("Variable") +
  geom_hline(yintercept = c(-threshold, threshold), col = "grey", linetype = 2)
  ylab("Estimates/CIs")


