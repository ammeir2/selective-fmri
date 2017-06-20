# Getting results ----------------------
simResults <- list()
slot <- 1
for(i in 200:400) {
  for(j in 1:60) {
    res <- NULL
    print(c(i, j))
    try(res <- readRDS(file = paste("fromcluster/simfmri may23 ", i," ", j,".RDS", sep ="")),
        silent = TRUE)
    if(!is.null(res)){
      simResults[[slot]] <- res
      slot <- slot + 1
    }
  }
}


# Preparing result list -------------------
library(ggplot2)
library(dplyr)
library(reshape2)
#simResults <- lapply(simResults, function(x) x[[1]])
resultList <- list()
len <- 1
for(i in 1:length(simResults)) {
  for(j in 1:length(simResults[[i]])) {
    nonEmpty <- sapply(simResults[[i]][[j]], length) > 0
    iterResult <- data.frame(matrix(nrow = sum(nonEmpty), ncol = 14))
    iterList <- simResults[[i]][[j]]
    row <- 1
    exp <- j * runif(1)
    for(k in which(nonEmpty)) {
      method <- iterList[[k]][[1]][[3]]
      if(method == 1) {
        method <- "onesided"
      } else if(method == 2) {
        method <- "selected"
      } else if (method == 3) {
        method <- "naive"
      }
      iterResult[row, ] <- c(experiment = exp,
                             cluster = row,
                             snr = as.numeric(iterList[[k]][[1]][[1]]),
                             spread = as.numeric(iterList[[k]][[1]][[2]]),
                             method = method,
                             rho = as.numeric(iterList[[k]][[1]][[4]]),
                             BHlevel = as.numeric(iterList[[k]][[1]][[5]]),
                             grp_size = as.numeric(iterList[[k]][[1]][[6]]),
                             size = as.numeric(iterList[[k]][[1]][[7]]),
                             true = as.numeric(iterList[[k]][[1]][[8]]),
                             pval = as.numeric(iterList[[k]][[2]][[3]]),
                             cover = as.numeric(iterList[[k]][[2]][[2]]),
                             naive = iterList[[k]][[3]][[2]],
                             cond = iterList[[k]][[3]][[1]])
      row <- row + 1
    }
    iterResult <- data.frame(iterResult)
    names(iterResult) <- c("experiment", "cluster", "snr", "spread", "method", "rho",
                           "pthreshold", "grp_size", "size", "true", "pval", "cover", "naive", "cond")
    resultList[[len]] <- iterResult
    len <- len + 1
  }
}
results <- data.table::rbindlist(resultList)
results <- as.data.frame(results)
simresults <- results

results <- simresults
results$cover <- as.numeric(as.character(results$cover))
results$pval <- as.numeric(as.character(results$pval))
results$true <- as.numeric(as.character(results$true))
results$snr <- as.numeric(as.character(results$snr))
results$experiment <- as.numeric(as.character(results$experiment))
results$size <- as.numeric(as.character(results$size))
results$cond <- as.numeric(as.character(results$cond))
results$naive <- as.numeric(as.character(results$naive))
results$pthreshold <- as.numeric(as.character(results$pthreshold))

# Cover Rate and Power ----------------
level <- 0.05
cover <- summarize(group_by(results, experiment, snr, method, rho, pthreshold),
                   cover = sum((cover > level) * size, na.rm = TRUE) / sum(size[!is.na(cover > level)]),
                   power = pval[which.max(true)] < level)
cover <- summarize(group_by(cover, snr, method, rho, pthreshold),
                   coversd = sd(cover) / sqrt(length(cover)),
                   cover = mean(cover),
                   powersd = sd(power, na.rm = TRUE) / sqrt(length(power)),
                   power = mean(power, na.rm = TRUE))

quantile <- qnorm(1 - 0.05 / nrow(cover))
ciwidth <- 0.1
forplot <- subset(cover, method == "selected")
scoverplot <- ggplot(subset(forplot, pthreshold %in% c(0.01, 0.001))) +
  geom_line(aes(x = snr, y = cover, col = method, linetype = method)) +
  geom_point(aes(x = snr, y = cover, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_hline(yintercept = 1 - level) + theme_bw() +
  geom_segment(aes(y = pmax(cover - coversd * quantile, 0), yend = pmin(cover + coversd * quantile, 1),
                   x = snr, xend = snr, col = method, linetype = method)) +
  geom_segment(aes(y = pmax(cover - coversd * quantile, 0), yend = pmax(cover - coversd * quantile, 0),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  geom_segment(aes(y = pmin(cover + coversd * quantile, 1), yend = pmin(cover + coversd * quantile, 1),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  ggtitle("Simulated Residuals") + theme(legend.position = "none") +
  scale_colour_brewer(palette = "Set1") + ylim(0.9, 1)
scoverplot

spowerplot <- ggplot(subset(forplot, pthreshold %in% c(0.01, 0.001))) +
  geom_line(aes(x = snr, y = power, col = method, linetype = method)) +
  geom_point(aes(x = snr, y = power, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_hline(yintercept = level) + theme_bw() +
  geom_segment(aes(y = pmax(power - powersd * quantile, 0), yend = pmin(power + powersd * quantile, 1),
                   x = snr, xend = snr, col = method, linetype = method)) +
  geom_segment(aes(y = pmax(power - powersd * quantile, 0), yend = pmax(power - powersd * quantile, 0),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  geom_segment(aes(y = pmin(power + powersd * quantile, 1), yend = pmin(power + powersd * quantile, 1),
                   x = snr - ciwidth, xend = snr + ciwidth, col = method)) +
  xlim(1 - ciwidth, 4 + ciwidth) + ylim(0, 0.8) +
  ggtitle("Simulated Residuals") + theme(legend.position = "none")
spowerplot

pdf("figures/ciplots.pdf",pagecentre=T, width=8, height=3 ,paper = "special")
gridExtra::grid.arrange(coverplot, powerplot, nrow = 1)
dev.off()


# RMSE and Bias -------------------
naive <- subset(results, method == "selected")
naive <- summarize(group_by(naive, snr, pthreshold, rho, experiment),
                   mse = weighted.mean(sqrt((naive - true)^2), size),
                   bias = weighted.mean(sign(naive) * (naive - true), size))
naive <- summarize(group_by(naive, snr, pthreshold, rho),
                   msesd = sd(mse) / sqrt(length(mse)), mse = mean(mse),
                   biassd = sd(bias) / sqrt(length(bias)), bias = mean(bias))
naive$method <- "naive"

cond <- subset(results, method == "selected")
cond <- summarize(group_by(cond, snr, pthreshold, rho, experiment),
                  mse = weighted.mean(sqrt((cond - true)^2), size),
                  bias = weighted.mean(sign(naive)*(cond - true), size))
cond <- summarize(group_by(cond, snr, pthreshold, rho),
                  msesd = sd(mse) / sqrt(length(mse)), mse = mean(mse),
                  biassd = sd(bias) / sqrt(length(bias)), bias = mean(bias))
cond$method <- "conditional"

forplot <- rbind(naive, cond)
smse <- ggplot(subset(forplot, pthreshold < 0.1), aes(x = snr, y = mse, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_line(aes(linetype = method)) + geom_point() +
  geom_segment(aes(y = mse - 2 * msesd, yend = mse + 2 * msesd,
                   x = snr, xend = snr, col = method), linetype = 2) +
  theme_bw() + ggtitle("Simulated Residuals") + ylab("RMSE") +
  theme(legend.position = "none")
smse

sbias <- ggplot(subset(forplot, pthreshold < 0.1), aes(x = snr, y = bias, col = method)) +
  facet_grid(rho ~ pthreshold, labeller = "label_both") +
  geom_line(aes(linetype = method)) + geom_point() +
  geom_segment(aes(y = bias - 2 * biassd, yend = bias + 2 * biassd,
                   x = snr, xend = snr, col = method), linetype = 2) +
  geom_hline(yintercept = 0) + theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Simulated Residuals")
sbias

pdf("figures/estplots.pdf",pagecentre=T, width=8, height=3 ,paper = "special")
gridExtra::grid.arrange(mse, bias, nrow = 1)
dev.off()











