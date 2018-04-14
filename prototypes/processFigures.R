library(magrittr)
library(reshape2)
library(ggplot2)

set.seed(3)
n <- 200
p <- 50
rho <- 0.4
sigma <- rho^(as.matrix(dist(1:p)))
X <- mvtnorm::rmvnorm(n, sigma = sigma) %>% scale()
sparsity <- 2
snr <- 0.3
true <- rep(0, p)
true[sample.int(p, sparsity)] <- rnorm(sparsity)
mu <- as.numeric(X %*% true)
ysig <- sqrt(var(mu) / snr)
y <- rnorm(n, sd = ysig) + mu
y <- y - mean(y)

# Raw Data
suffStat <- as.numeric(t(X) %*% y / n)
X <- X[, order(abs(suffStat), decreasing = TRUE)]
suffStat <- suffStat[order(abs(suffStat), decreasing = TRUE)]
suffDat <- data.frame(variable = 1:p, suffStat = suffStat)
ggplot(suffDat) +
  geom_point(aes(x = variable, y = suffStat)) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  ylab("Sufficient Statistics") +
  xlab("Variable") +
  ggtitle("Data") +
  ylim(-0.75, 0.3)
ggsave("vartex/demo_raw.pdf", width = 4, height = 3)

# W selected
selected <- 1:p %in% order(abs(suffStat), decreasing = TRUE)[1:10]
suffDat$selected <- selected
ggplot(suffDat) +
  geom_point(aes(x = variable, y = suffStat, col = selected)) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  ylab("Sufficient Statistics") +
  xlab("Variable") +
  ggtitle("Data -> Selection") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 10.5, linetype = 2) +
  ylim(-0.75, 0.3)
ggsave("vartex/demo_selection.pdf", width = 4, height = 3)

select <- (1:p) %in% order(abs(suffStat), decreasing = TRUE)[1:10]
threshold <- (min(abs(suffStat[select])) + max(abs(suffStat)[!select])) / 2 * n
fit <- approxConditionalMLE(X, y, ysig, threshold,
                            thresholdLevel = 0.01, cilevel = 0.05,
                            varCI = FALSE, bootSamples = 2000,
                            verbose = TRUE, computeMLE = TRUE)
fit$sampMean


# Thresholding
suffDat$sampMean <- fit$sampMean / n
ggplot(suffDat) +
  geom_point(aes(x = variable, y = sampMean, col = selected)) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  ylab("Sufficient Statistics") +
  xlab("Variable") +
  ggtitle("Data -> Selection -> Thresholding") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 10.5, linetype = 2) +
  ylim(-0.75, 0.3)
ggsave("vartex/demo_thresholding.pdf", width = 4, height = 3)

# Thresholded Coefficients
suffDat$sampMean <- fit$sampMean / n
suffDat$sampCoef <- 0
suffDat$sampCoef[fit$sampCoef != 0] <- fit$sampCoef
ggplot(subset(suffDat, select)) +
  geom_point(aes(x = variable, y = sampCoef), col = "cyan3", size = 2) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  ylab("Sufficient Statistics") +
  xlab("Variable") +
  ggtitle("Thresholded Coefficients") +
  theme(legend.position = "none")
ggsave("vartex/demo_thresCoef.pdf", width = 4, height = 3)

# Thresholded Distribution
bootdist <- fit$suffSamp / n
bootdist <- melt(bootdist)
names(bootdist) <- c("replication", "variable", "sample")
bootdist$selected <- FALSE
bootdist$selected <- bootdist$variable %in% which(select)
ggplot(bootdist) +
  geom_violin(aes(x = factor(variable), y = sample, col = selected)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Variable") + ylab("Bootstrapped Sufficient Statistics") +
  ggtitle("Bootstrap Distribution")
ggsave("vartex/demo_suffBoot.pdf", width = 4, height = 3)


# Coefficient Sample
coefDist <- fit$naiveBoot
coefDist <- melt(coefDist)
names(coefDist) <- c("replication", "variable", "sample")
ggplot(coefDist) +
  geom_violin(aes(x = factor(variable), y = sample), col = "cyan3") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Variable") + ylab("Bootstrapped Coefficients") +
  ggtitle("Bootstrapped Coefficients") +
  ylim(-1.1, .7)
ggsave("vartex/demo_bootCoef.pdf", width = 4, height = 3)

# Residuals
coefDist <- fit$naiveBoot
for(j in 1:ncol(coefDist)) {
  coefDist[, j] <- coefDist[, j] - fit$sampCoef[j]
}
coefDist <- melt(coefDist)
names(coefDist) <- c("replication", "variable", "sample")
ggplot(coefDist) +
  geom_violin(aes(x = factor(variable), y = sample), col = "cyan3") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Variable") + ylab("Bootstrapped Coefficients") +
  ggtitle("Bootstrapped Coefficients -> Residuals")  +
  ylim(-1.1, .7)
ggsave("vartex/demo_residBoot.pdf", width = 4, height = 3)

# Coefficients -------
naivefit <- lm(y ~ X[, select] - 1)
sds <- sqrt(diag(vcov(naivefit)))
naiveCI <- matrix(nrow = sum(select), ncol = 2)
naiveCI[, 1] <- coef(naivefit) - 2 * sds
naiveCI[, 2] <- coef(naivefit) + 2 * sds
trueCoef <- coef(lm(mu ~ X[, select] - 1))
plotEstimates(estimates = list(naive = fit$naive, mle = fit$mle, true = trueCoef),
              CIs = list(naive = naiveCI, bootstrap = fit$naiveBootCI, poly = fit$polyCI))
ggsave(file = "vartex/demo_ci.pdf", width = 7, height = 3)

plotEstimates(estimates = list(naive = fit$naive, mle = fit$mle, true = trueCoef),
              CIs = list(naive = naiveCI, bootstrap = fit$naiveBootCI, poly = fit$polyCI)) +
  ylim(-1, 0.6)
ggsave(file = "vartex/demo_ci_nopoly.pdf", width = 7, height = 3)






