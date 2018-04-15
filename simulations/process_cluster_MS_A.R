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
# filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_contrastZ_A_*"))
# filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_univZ_A_*"))
# filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_univConstZ_F_*"))
filenames <- as.list(dir(path = 'simulations/results', pattern="variationalSim_oneStep_A_*"))
filenames <- lapply(filenames, function(x) paste0('simulations/results/', x))
results <- lapply(filenames, function(x) try(readRDS(x)))
results <- do.call("c", results)
results <- do.call("c", results)

library(dplyr)
library(reshape2)
library(ggplot2)
# Cover ---
computeCover <- function(x) {
  cis <- x$cis
  naive <- x$cis$naive
  # cis[2:3] <- lapply(cis[2:3], function(x) cbind(pmin(x[, 1], naive[, 1]), pmax(x[, 2], naive[, 2])))
  x$cis <- cis
  x$cis$mleCI <- cbind(pmin(x$cis$mleCI[, 1], naive[, 1]), pmax(x$cis$mleCI[, 2], naive[, 2]))

  cover <- sapply(x$cis[-6], getCover, x$estimate[, 5])
  result <- c(x$config, cover)
  # print(length(result))
  return(result)
}
coverdat <- lapply(results, function(x) try(computeCover(x)))
table(sapply(coverdat, length))
coverdat <- coverdat[sapply(coverdat, length) == 13]
coverdat <- do.call("rbind", coverdat)
coverdat <- data.frame(coverdat)

coverdat$reps <- NULL
coverdat <- melt(coverdat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(coverdat)[8:9] <- c("type", "cover")
coverdat <- summarize(group_by(coverdat, n, p, snr, sparsity, covtype, nselect, type, rho),
                      sd = sd(cover, na.rm = TRUE) / sqrt(length(cover)) ,
                      cover = mean(cover, na.rm = TRUE))
quant <- qnorm(1 - 0.05 / 2)
coverdat$type <- as.character(coverdat$type)
# coverdat$type[coverdat$type == "naiveBoot"] <- "bootstrap"
# coverdat <- subset(coverdat, covtype == 2 & snr != 0.01 & snr <= 2 & type %in% c("naive", "bootstrap", "poly"))
# coverdat$type <- factor(coverdat$type, levels = c("naive", "bootstrap", "poly"))
coverdat$rholab <- "Independent Design"
coverdat$rholab[coverdat$rho == 0.7] <- "AR Design"
coverdat$rholab <- factor(coverdat$rholab, levels = c("Independent Design", "AR Design"))
coverdat$sparselab <- paste("Sparsity: ", coverdat$sparsity, sep = "")
ggplot(coverdat) +
  geom_line(aes(x = log2(snr), y = cover, col = type, linetype = type)) +
  facet_grid(rholab ~ sparselab) +
  theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = pmin(cover + quant * sd, 1), yend = pmax(cover - quant*sd, 0))) +
  geom_hline(yintercept = 0.95) +
  ylab("Coverage Rate") + xlab("log2(snr)")
# ggsave("vartex/coversim.pdf", height = 5, width = 7)

# Plot for Northshore -------
ggplot(subset(coverdat, sparsity == 2 & rholab == "Independent Design" & type == "naive")) +
  geom_line(aes(x = log2(snr), y = cover, col = type, linetype = type)) +
  theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = pmin(cover + quant * sd, 1), yend = pmax(cover - quant*sd, 0))) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  ylab("Coverage Rate") + xlab("log2(snr)")
ggsave(file = "/Users/amitmeir/Documents/academics/northshore presentation/naiveCoverMS.pdf",
       width = 5, height = 3)

coverdat$type <- as.character(coverdat$type)
coverdat$type[coverdat$type == "bootstrap"] <- "adjusted"
coverdat$type <- factor(coverdat$type, levels = c("naive", "adjusted", "poly"))
ggplot(subset(coverdat, sparsity == 2 & rholab == "Independent Design" & type != "poly")) +
  geom_line(aes(x = log2(snr), y = cover, col = type, linetype = type)) +
  theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = pmin(cover + quant * sd, 1), yend = pmax(cover - quant*sd, 0))) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  ylab("Coverage Rate") + xlab("log2(snr)")
ggsave(file = "/Users/amitmeir/Documents/academics/northshore presentation/adjCoverMS.pdf",
       width = 5, height = 3)


# Power ---
computePower <- function(x) {
  cis <- x$cis
  naive <- x$cis$naive
  # cis[2:3] <- lapply(cis[2:3], function(x) cbind(pmin(x[, 1], naive[, 1]), pmax(x[, 2], naive[, 2])))
  x$cis <- cis
  x$cis$mleCI <- cbind(pmin(x$cis$mleCI[, 1], naive[, 1]), pmax(x$cis$mleCI[, 2], naive[, 2]))

  nonzero <- x$estimate$true != 0
  if(all(!nonzero)) {
    return(c(x$config, sapply(cis, function(x) NA)))
  }

  cbind(cis$varBoot, cis$hybrid, cis$poly, x$estimate$true)
  power <- sapply(cis, function(ci, nonzero) mean(sign(ci[nonzero, 1]) == sign(ci[nonzero, 2])), nonzero)
  result <- c(x$config, power)
  return(result)
}
powerdat <- lapply(results, function(x) try(computePower(x)))
table(sapply(powerdat, length))
powerdat <- powerdat[which(sapply(powerdat, length) == 14)]
powerdat <- do.call("rbind", powerdat)
powerdat <- data.frame(powerdat)
powerdat$reps <- NULL
powerdat <- melt(powerdat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(powerdat)[8:9] <- c("type", "power")
powerdat <- summarize(group_by(powerdat, n, p, snr, sparsity, covtype, nselect, type, rho),
                      sd = sd(power, na.rm = TRUE) / sqrt(length(power)) ,
                      power = mean(power, na.rm = TRUE))
powerdat$type <- as.character(powerdat$type)
powerdat$type[powerdat$type == "naiveBoot"] <- "bootstrap"
powerdat$type <- factor(powerdat$type, levels = c("naive", "bootstrap", "poly"))
powerdat$sparselab <- paste("Sparsity:", powerdat$sparsity)
powerdat$rholab <- "Independent Design"
powerdat$rholab[powerdat$rho == 0.7] <- "AR Design"
powerdat$rholab <- factor(powerdat$rholab, levels = c("Independent Design", "AR Design"))
ggplot(subset(powerdat, snr != 0.01 & snr <= 2 &
                type %in% c("naive", "poly", "bootstrap"))) +
  geom_line(aes(x = log2(snr), y = power, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = power, col = type, shape = type)) +
  facet_grid(rholab ~ sparselab, scales = "free_y") +
  theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = pmin(power + quant * sd, 1), yend = pmax(power - quant*sd, 0))) +
  ylim(0, 1) + ylab("Power")
# ggsave("vartex/powersim.pdf", height = 5, width = 7)

# Size -----
computeSize <- function(x) {
  cis <- x$cis
  naive <- x$cis$naive
  # cis[2:3] <- lapply(cis[2:3], function(x) cbind(pmin(x[, 1], naive[, 1]), pmax(x[, 2], naive[, 2])))
  x$cis <- cis
  x$cis$mleCI <- cbind(pmin(x$cis$mleCI[, 1], naive[, 1]), pmax(x$cis$mleCI[, 2], naive[, 2]))

  size <- sapply(x$cis, function(x) median(x[, 2] - x[, 1]))
  size <- log2(size[-1] / size[1])
  result <- c(x$config, size)
  return(result)
}
sizedat <- lapply(results, function(x) try(computeSize(x)))
table(sapply(sizedat, length))
sizedat <- sizedat[which(sapply(sizedat, length) == 13)]
sizedat <- do.call("rbind", sizedat)
sizedat <- data.frame(sizedat)
sizedat$reps <- NULL
sizedat <- melt(sizedat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(sizedat)[8:9] <- c("type", "relsize")
sizedat <- summarize(group_by(sizedat, n, p, snr, sparsity, covtype, nselect, type, rho),
                     sd = sd(relsize, na.rm = TRUE) / sqrt(length(relsize)) ,
                     relsize = mean(relsize, na.rm = TRUE))
sizedat$type <- as.character(sizedat$type)
sizedat$type[sizedat$type == "naiveBoot"] <- "bootstrap"
sizedat$type <- factor(sizedat$type, levels = c("bootstrap", "poly"))
sizedat$sparselab <- paste("Sparsity:", sizedat$sparsity)
sizedat$rholab <- "Independent Design"
sizedat$rholab[sizedat$rho == 0.7] <- "AR Design"
sizedat$rholab <- factor(sizedat$rholab, levels = c("Independent Design", "AR Design"))
ggplot(subset(sizedat, snr != 0.01 & snr <= 2 &
                type %in% c("bootstrap", "poly", "naive"))) +
  geom_line(aes(x = log2(snr), y = relsize, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = relsize, col = type, shape = type)) +
  facet_grid(rholab ~ sparselab) + theme_bw() +
  geom_segment(aes(x = log2(snr), xend = log2(snr), col = type, linetype = type,
                   y = relsize + quant * sd, yend = relsize - quant*sd)) +
  geom_hline(yintercept = 0) +
  ylab("log2(CI size / naive CI size)")
# ggsave("vartex/sizesim.pdf", height = 5, width = 7)

# Rel MSE -------
computeRelMSE <- function(x) {
  true <- x$estimate[, 6]
  mse <- apply(x$estimate[, -c(5:6)], 2, function(x) sqrt(mean((x - true)^2)))
  relmse <- log2(mse[1:3] / mse[4])
  result <- c(x$config, relmse)
  return(result)
}
relmsedat <- lapply(results, function(x) try(computeRelMSE(x)))
table(sapply(relmsedat, length))
relmsedat <- relmsedat[sapply(relmsedat, length) == 11]
relmsedat <- do.call("rbind", relmsedat)
relmsedat <- data.frame(relmsedat)
relmsedat$reps <- NULL
relmsedat <- melt(relmsedat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
names(relmsedat)[8:9] <- c("type", "rmse")
relmsedat <- summarize(group_by(relmsedat, n, p, snr, sparsity, covtype, nselect, type, rho),
                       sd = sd(rmse, na.rm = TRUE) / sqrt(length(rmse)) ,
                       rmse = mean(rmse, na.rm = TRUE))
relmsedat$rholab <- "Independent Design"
relmsedat$rholab[relmsedat$rho == 0.7] <- "AR Design"
relmsedat$rholab <- factor(relmsedat$rholab, levels = c("Independent Design", "AR Design"))
relmsedat$sparselab <- "Sparsity = 2"
relmsedat$sparselab[relmsedat$sparsity == 8] <- "Sparsity = 8"
ggplot(subset(relmsedat, covtype == 2 & type != "var")) +
  geom_line(aes(x = log2(snr), y = rmse, col = type, linetype = type)) +
  geom_point(aes(x = log2(snr), y = rmse, col = type, shape = type)) +
  facet_grid(rholab ~ sparselab) + theme_bw() +
  geom_hline(yintercept = 0) +
  geom_segment(aes(x = log2(snr), xend = log2(snr),
                   y = rmse - quant*sd, yend = rmse + quant*sd,
                   linetype = type, col = type)) +
  ylab("Relative MSE") + xlab("log2(snr)")
ggsave(filen = "vartex/relMSE.pdf", width = 7, height = 5)

# MSE ----
# computeMSE <- function(x) {
#   true <- x$estimate[, 6]
#   mse <- apply(x$estimate[, -c(5:6)], 2, function(x) sqrt(mean((x - true)^2)))
#   result <- c(x$config, mse)
#   return(result)
# }
# msedat <- t(sapply(results, computeMSE))
# msedat <- data.frame(msedat)
# msedat$reps <- NULL
# msedat <- melt(msedat, id = c("n", "p", "snr", "sparsity", "covtype", "nselect", "rho"))
# names(msedat)[8:9] <- c("type", "rmse")
# msedat <- summarize(group_by(msedat, n, p, snr, sparsity, covtype, nselect, type, rho),
#                     sd = sd(rmse, na.rm = TRUE) / sqrt(length(rmse)) ,
#                     rmse = mean(rmse, na.rm = TRUE))
# ggplot(subset(msedat, covtype == 2)) +
#   geom_point(aes(x = log2(snr), y = rmse, col = type, shape = type)) +
#   geom_line(aes(x = log2(snr), y = rmse, col = type, linetype = type)) +
#   facet_grid(rho ~ sparsity, labeller = "label_both") + theme_bw()
#
