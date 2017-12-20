getCover <- function(ci, truth) {
  cover <- 0
  for(i in 1:length(truth)) {
    if(ci[i, 1] < truth[i] & ci[i, 2] > truth[i]) {
      cover <- cover + 1
    } else {
      #cat(i, " ")
    }
  }
  return(cover / length(truth))
}

# filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_H_*"))
filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_contrastZ_A_*"))
filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_univZ_A_*"))
filenames <- lapply(filenames, function(x) paste0('simulations/results/', x))[-c(3, 4)]
results <- lapply(filenames, function(x) try(readRDS(x)))
results <- do.call("c", results)
results <- do.call("c", results)

library(dplyr)
library(reshape2)
library(ggplot2)
# Cover ---
computeCover <- function(x) {
  cis <- x$cis
  # naive <- x$cis$naive
  # cis[2:3] <- lapply(cis[2:3], function(x) cbind(pmin(x[, 1], naive[, 1]), pmax(x[, 2], naive[, 2])))
  x$cis <- cis

  cover <- sapply(x$cis[-6], getCover, x$estimate[, 5])
  result <- c(x$config, cover)
  return(result)
}
coverdat <- t(sapply(results, computeCover))
coverdat <- data.frame(coverdat)
coverdat$reps <- NULL
coverdat <- melt(coverdat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(coverdat)[8:9] <- c("type", "cover")
coverdat <- summarize(group_by(coverdat, n, p, snr, sparsity, covtype, nselect, type, rho),
                      sd = sd(cover, na.rm = TRUE) / sqrt(length(cover)) ,
                      cover = mean(cover, na.rm = TRUE))
ggplot(subset(coverdat, covtype == 2 & snr != 0.01 & snr <= 1)) +
  geom_line(aes(x = log2(snr), y = cover, col = type, linetype = type)) +
  facet_grid(rho ~ sparsity, labeller = "label_both") + theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type,
                   y = cover + 2 * sd, yend = cover - 2*sd)) +
  geom_hline(yintercept = 0.95)

# Power ---
computePower <- function(x) {
  cis <- x$cis
  # naive <- x$cis$naive
  # cis[2:3] <- lapply(cis[2:3], function(x) cbind(pmin(x[, 1], naive[, 1]), pmax(x[, 2], naive[, 2])))
  x$cis <- cis

  nonzero <- x$estimate$true != 0
  if(all(!nonzero)) {
    return(c(x$config, sapply(cis, function(x) NA)))
  }

  cbind(cis$varBoot, cis$hybrid, cis$poly, x$estimate$true)
  power <- sapply(cis, function(ci, nonzero) mean(sign(ci[nonzero, 1]) == sign(ci[nonzero, 2])), nonzero)
  result <- c(x$config, power)
  return(result)
}
powerdat <- lapply(results, computePower)
# powerdat <- powerdat[-which(sapply(powerdat, length) == 14)]
powerdat <- do.call("rbind", powerdat)
powerdat <- data.frame(powerdat)
powerdat$reps <- NULL
powerdat <- melt(powerdat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(powerdat)[8:9] <- c("type", "power")
powerdat <- summarize(group_by(powerdat, n, p, snr, sparsity, covtype, nselect, type, rho),
                      sd = sd(power, na.rm = TRUE) / sqrt(length(power)) ,
                      power = mean(power, na.rm = TRUE))
ggplot(subset(powerdat, snr != 0.01 & snr <= 2)) +
  geom_line(aes(x = log2(snr), y = power, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = power, col = type, shape = type)) +
  facet_grid(rho ~ sparsity, labeller = "label_both", scales = "free_y") +
  theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = power + 2 * sd, yend = power - 2*sd)) +
  ylim(0, 1)


# Size -----
computeSize <- function(x) {
  cis <- x$cis
  # naive <- x$cis$naive
  # cis[2:3] <- lapply(cis[2:3], function(x) cbind(pmin(x[, 1], naive[, 1]), pmax(x[, 2], naive[, 2])))
  x$cis <- cis

  size <- sapply(x$cis, function(x) median(x[, 2] - x[, 1]))
  size <- log2(size[-1] / size[1])
  result <- c(x$config, size)
  return(result)
}
sizedat <- lapply(results, computeSize)
# sizedat <- sizedat[-which(sapply(sizedat, length) == 13)]
sizedat <- do.call("rbind", sizedat)
sizedat <- data.frame(sizedat)
sizedat$reps <- NULL
sizedat <- melt(sizedat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(sizedat)[8:9] <- c("type", "relsize")
sizedat <- summarize(group_by(sizedat, n, p, snr, sparsity, covtype, nselect, type, rho),
                     sd = sd(relsize, na.rm = TRUE) / sqrt(length(relsize)) ,
                     relsize = mean(relsize, na.rm = TRUE))
ggplot(subset(sizedat, snr != 0.01 & snr <= 2 & type != "ind")) +
  geom_line(aes(x = log2(snr), y = relsize, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = relsize, col = type, shape = type)) +
  facet_grid(rho ~ sparsity, labeller = "label_both") + theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = relsize + 2 * sd, yend = relsize - 2*sd)) +
  geom_hline(yintercept = 0)

# MSE ----
computeMSE <- function(x) {
  true <- x$estimate[, 6]
  mse <- apply(x$estimate[, -c(5:6)], 2, function(x) sqrt(mean((x - true)^2)))
  result <- c(x$config, mse)
  return(result)
}
msedat <- t(sapply(results, computeMSE))
msedat <- data.frame(msedat)
msedat$reps <- NULL
msedat <- melt(msedat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(msedat)[8:9] <- c("type", "rmse")
msedat <- summarize(group_by(msedat, n, p, snr, sparsity, covtype, nselect, type, rho),
                    sd = sd(rmse, na.rm = TRUE) / sqrt(length(rmse)) ,
                    rmse = mean(rmse, na.rm = TRUE))
ggplot(subset(msedat, covtype == 2)) +
  geom_point(aes(x = log2(snr), y = rmse, col = type, shape = type)) +
  geom_line(aes(x = log2(snr), y = rmse, col = type, linetype = type)) +
  facet_grid(rho ~ sparsity, labeller = "label_both") + theme_bw()

# Rel MSE -------
computeRelMSE <- function(x) {
  true <- x$estimate[, 6]
  mse <- apply(x$estimate[, -c(5:6)], 2, function(x) sqrt(mean((x - true)^2)))
  relmse <- log2(mse[1:3] / mse[4])
  result <- c(x$config, relmse)
  return(result)
}
relmsedat <- t(sapply(results, computeRelMSE))
relmsedat <- data.frame(relmsedat)
relmsedat$reps <- NULL
relmsedat <- melt(relmsedat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(relmsedat)[8:9] <- c("type", "rmse")
relmsedat <- summarize(group_by(relmsedat, n, p, snr, sparsity, covtype, nselect, type, rho),
                    sd = sd(rmse, na.rm = TRUE) / sqrt(length(rmse)) ,
                    rmse = mean(rmse, na.rm = TRUE))
ggplot(subset(relmsedat, covtype == 2 & snr != 0.01)) +
  geom_line(aes(x = log2(snr), y = rmse, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = rmse, col = type, shape = type)) +
  facet_grid(rho ~ sparsity, labeller = "label_both") + theme_bw() +
  geom_hline(yintercept = 0) +
  geom_segment(aes(x = log2(snr), xend = log2(snr), y = rmse - 2*sd, yend = rmse + 2*sd,
                   linetype = type, col = type)) +
  ylab("Relative MSE") + xlab("log2(sample size)")


