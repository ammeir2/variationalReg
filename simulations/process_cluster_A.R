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

filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_A*"))
filenames <- lapply(filenames, function(x) paste0('simulations/results/', x))[-c(3, 4)]
results <- lapply(filenames, function(x) readRDS(x))
results <- do.call("c", results)
results <- do.call("c", results)

library(dplyr)
library(reshape2)
library(ggplot2)
# MSE ----
computeMSE <- function(x) {
  true <- x$estimate[, 5]
  mse <- apply(x$estimate[, -c(4:5)], 2, function(x) sqrt(mean((x - true)^2)))
  result <- c(x$config, mse)
  return(result)
}
msedat <- t(sapply(results, computeMSE))
msedat <- data.frame(msedat)
msedat$reps <- NULL
msedat <- melt(msedat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect"))
names(msedat)[7:8] <- c("type", "rmse")
msedat <- summarize(group_by(msedat, n, p, snr, sparsity, covtype, nselect, type),
                    sd = sd(rmse, na.rm = TRUE) / sqrt(length(rmse)) ,
                    rmse = mean(rmse, na.rm = TRUE))
ggplot(subset(msedat, covtype == 3)) + geom_line(aes(x = log2(n), y = rmse, col = type, linetype = type)) +
  facet_grid(sparsity ~ snr, labeller = "label_both") + theme_bw()

# Cover ---
computeCover <- function(x) {
  cover <- sapply(x$cis, getCover, x$estimate[, 5])
  result <- c(x$config, cover)
  return(result)
}
coverdat <- t(sapply(results, computeCover))
coverdat <- data.frame(coverdat)
coverdat$reps <- NULL
coverdat <- melt(coverdat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect"))
names(coverdat)[7:8] <- c("type", "cover")
coverdat <- summarize(group_by(coverdat, n, p, snr, sparsity, covtype, nselect, type),
                    sd = sd(cover, na.rm = TRUE) / sqrt(length(cover)) ,
                    cover = mean(cover, na.rm = TRUE))
ggplot(subset(coverdat, covtype == 1)) +
  geom_line(aes(x = log2(n), y = cover, col = type, linetype = type)) +
  facet_grid(sparsity ~ snr, labeller = "label_both") + theme_bw() +
  geom_segment(aes(x = log2(n), xend = log2(n), col = type,
                   y = cover + 2 * sd, yend = cover - 2*sd)) +
  geom_hline(yintercept = 0.95)

