sampleTrunc <- function(n, mu = 0, sd = 1) {
  nsamp <- 1
  samp <- numeric(n)
  threshold <- 1.5
  tries <- 0
  while(nsamp <= n) {
    tries <- tries + 1
    y <- rnorm(1, mu, sd)
    if(abs(y) > threshold) {
      samp[nsamp] <- y
      nsamp <- nsamp + 1
    }
  }
  return(samp)
}

library(dplyr)
library(ggplot2)
library(reshape2)

# Naive Bootstrap -------------------
n <- 10^3
samp <- sampleTrunc(n)
bootDist <- sapply(samp, function(y) sampleTrunc(1000, mu = y, sd = 1) - y)
trueDist <- sampleTrunc(10^4, mu = 0, sd = 1)
residquants <- apply(bootDist, 2, quantile, c(0.025, 0.975))
cis <- matrix(nrow = n, ncol = 3)
for(i in 1:nrow(cis)) {
  cis[i, 1] <- samp[i]
  cis[i, 2] <- samp[i] - residquants[2, i]
  cis[i, 3] <- samp[i] - residquants[1, i]
}
cis <- data.frame(cis)
names(cis) <- c("estimate", "lci", "uci")
cis$cover <- apply(cis, 1, function(x) x[2] < 0 & x[3] > 0)
mean(cis[, 3] - cis[, 2])
ggplot(cis) +
  geom_segment(aes(x = estimate, xend = estimate, y = lci, yend = uci, col = cover), alpha = 0.5) +
  theme_bw() + geom_hline(yintercept = 0) +
  xlab("Observed Y") + ylab("Confidence Intervals") +
  ggtitle("Naive Bootstrap CIs") +
  theme(legend.position = "none")
ggsave(file = "vartex/naiveBootCIcover.pdf", width = 3, height = 3)


bootDist <- melt(bootDist[, 1:3])
names(bootDist) <- c("ind", "muhat", "resid")
bootDist$muhat <- round(y[bootDist$muhat], 2)
trueDist <- data.frame(ind = 1, muhat = 0, resid = trueDist)
forplot <- rbind(bootDist, trueDist)
muvals <- unique(forplot$muhat)
forplot$muhat <- factor(forplot$muhat, levels = c(0, subset(muvals, muvals != 0)))
forplot$mu <- forplot$muhat
ggplot(forplot) + geom_density(aes(x = resid, col = mu)) +
  theme_bw() + facet_wrap(~ mu, labeller = "label_both") +
  ylab("Density of Residual Distribution") + xlab("") +
  theme(legend.position = "none") + ggtitle("Bootstrap Distribution")
ggsave(file = "vartex/naiveBootDist.pdf", width = 4, height = 3)


# Conditional Bootstrap -----------
n <- 1000
samp <- sampleTrunc(n)
mutild <- samp * (1 - pnorm(abs(samp)) < 0.0025)
bootDist <- sapply(mutild, function(est) sampleTrunc(1000, mu = est, sd = 1) - est)
trueDist <- sampleTrunc(10^4, mu = 0, sd = 1)
residquants <- apply(bootDist, 2, quantile, c(0.025, 0.975))
cis <- matrix(nrow = n, ncol = 3)
for(i in 1:nrow(cis)) {
  cis[i, 1] <- samp[i]
  cis[i, 2] <- samp[i] - residquants[2, i]
  cis[i, 3] <- samp[i] - residquants[1, i]
}
cis <- data.frame(cis)
names(cis) <- c("estimate", "lci", "uci")
cis$cover <- apply(cis, 1, function(x) x[2] < 0 & x[3] > 0)
mean(cis[, 3] - cis[, 2])
ggplot(cis) +
  geom_segment(aes(x = estimate, xend = estimate, y = lci, yend = uci, col = cover), alpha = 0.5) +
  theme_bw() + geom_hline(yintercept = 0) +
  xlab("Observed Y") + ylab("Confidence Intervals") +
  ggtitle("Conditional Bootstrap CIs")
ggsave(file = "vartex/condUnivBootCIcover.pdf", width = 6, height = 3)
mean(cis$cover)


bootDist <- melt(bootDist[, 1:3])
names(bootDist) <- c("ind", "muhat", "resid")
bootDist$muhat <- round(y[bootDist$muhat], 2)
trueDist <- data.frame(ind = 1, muhat = 0, resid = trueDist)
forplot <- rbind(bootDist, trueDist)
muvals <- unique(forplot$muhat)
forplot$muhat <- factor(forplot$muhat, levels = c(0, subset(muvals, muvals != 0)))
forplot$mu <- forplot$muhat
ggplot(forplot) + geom_density(aes(x = resid, col = mu)) +
  theme_bw() + facet_wrap(~ mu, labeller = "label_both") +
  ylab("Density of Residual Distribution") + xlab("") +
  theme(legend.position = "none") + ggtitle("Bootstrap Distribution")
ggsave(file = "vartex/condUnivBootDist.pdf", width = 4, height = 3)








