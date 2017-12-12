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

filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSimAgg_B_*"))
filenames <- lapply(filenames, function(x) paste0('simulations/results/', x))[-c(3, 4)]
results <- lapply(filenames, function(x) readRDS(x))
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

  cover <- sapply(x$cis, getCover, x$estimates$true)
  result <- c(x$config, cover)
  return(result)
}
coverdat <- t(sapply(results, computeCover))
coverdat <- data.frame(coverdat)
coverdat$reps <- NULL
coverdat <- melt(coverdat, id = c("n", "p", "snr", "sparsity", "covtype", "rho", "pthreshold"))
names(coverdat)[8:9] <- c("type", "cover")
coverdat <- summarize(group_by(coverdat, n, p, snr, sparsity, covtype, type, rho, pthreshold),
                    sd = sd(cover, na.rm = TRUE) / sqrt(length(cover)) ,
                    cover = mean(cover, na.rm = TRUE))
ggplot(subset(coverdat)) +
  geom_line(aes(x = log2(snr), y = cover, col = type, linetype = type)) +
  facet_grid(sparsity ~ pthreshold, labeller = "label_both") + theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type,
                   y = cover + 2 * sd, yend = cover - 2*sd)) +
  geom_hline(yintercept = 0.95)

# Power ---
computePower <- function(x) {
  cis <- x$cis
  # naive <- x$cis$naive
  # cis[2:3] <- lapply(cis[2:3], function(x) cbind(pmin(x[, 1], naive[, 1]), pmax(x[, 2], naive[, 2])))
  x$cis <- cis

  nonzero <- x$estimates$true != 0
  if(all(!nonzero)) {
    return(c(x$config, sapply(cis, function(x) NA)))
  }

  cbind(cis$varBoot, cis$hybrid, cis$poly, x$estimate$true)
  power <- sapply(cis, function(ci, nonzero) mean(sign(ci[nonzero, 1]) == sign(ci[nonzero, 2])), nonzero)
  result <- c(x$config, power)
  return(result)
}
powerdat <- t(sapply(results, computePower))
powerdat <- data.frame(powerdat)
powerdat$reps <- NULL
powerdat <- melt(powerdat, id = c("n", "p", "snr", "sparsity", "covtype", "pthreshold", "rho"))
names(powerdat)[8:9] <- c("type", "power")
powerdat <- summarize(group_by(powerdat, n, p, snr, sparsity, covtype, pthreshold, type, rho),
                     sd = sd(power, na.rm = TRUE) / sqrt(length(power)) ,
                     power = mean(power, na.rm = TRUE))
ggplot(subset(powerdat)) +
  geom_line(aes(x = log2(snr), y = power, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = power, col = type, shape = type)) +
  facet_grid(sparsity ~ pthreshold, labeller = "label_both", scales = "free_y") +
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

  size <- sapply(x$cis, function(x) mean(x[, 2] - x[, 1]))
  size <- log2(size[-1] / size[1])
  result <- c(x$config, size)
  return(result)
}
sizedat <- t(sapply(results, computeSize))
sizedat <- data.frame(sizedat)
sizedat$reps <- NULL
sizedat <- melt(sizedat, id = c("n", "p", "snr", "sparsity", "covtype", "rho", "pthreshold"))
names(sizedat)[8:9] <- c("type", "relsize")
sizedat <- summarize(group_by(sizedat, n, p, snr, sparsity, covtype, pthreshold, type, rho),
                     sd = sd(relsize, na.rm = TRUE) / sqrt(length(relsize)) ,
                     relsize = mean(relsize, na.rm = TRUE))
ggplot(subset(sizedat)) +
  geom_line(aes(x = log2(snr), y = relsize, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = relsize, col = type, shape = type)) +
  facet_grid(sparsity ~ pthreshold, labeller = "label_both") + theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = relsize + 2 * sd, yend = relsize - 2*sd)) +
  geom_hline(yintercept = 0)

