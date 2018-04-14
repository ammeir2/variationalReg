library(magrittr)
library(reshape2)
library(ggplot2)

set.seed(4)
# set.seed(5)
n <- 200
p <- 50
X <- rnorm(n * p) %>% matrix(ncol = p) %>% scale()
sparsity <- 3
snr <- 0.2
true <- rep(0, p)
true[sample.int(p, sparsity)] <- rnorm(sparsity)
mu <- as.numeric(X %*% true)
ysig <- sqrt(var(mu) / snr)
y <- rnorm(n, sd = ysig) + mu
y <- y - mean(y)

suffStat <- t(X) %*% y %>% as.numeric()
select <- (1:p) %in% order(abs(suffStat), decreasing = TRUE)[1:10]
threshold <- (min(abs(suffStat[select])) + max(abs(suffStat)[!select])) / 2
fit <- approxConditionalMLE(X, y, ysig, threshold,
                            thresholdLevel = 0.01, cilevel = 0.05,
                            varCI = FALSE, bootSamples = 2000,
                            verbose = TRUE, computeMLE = TRUE)

gradPath <- fit$sampPath[, select]
for(i in 1:nrow(gradPath)) {
  gradPath[i, ] <- (suffStat[select] - gradPath[i, ]) #/ i^0.6
}
gradPath <- apply(gradPath, 2, function(x) cumsum(x) / 1:length(x))

gradPath <- melt(gradPath)
names(gradPath) <- c("Iteration", "Variable", "Gradient")
gradPath <- cbind(gradPath, Value = "Gradient Cumulative Mean")

coefPath <- fit$coefPath
coefPath <- melt(coefPath)
names(coefPath) <- c("Iteration", "Variable", "Gradient")
coefPath <- cbind(coefPath, Value = "Coefficient Estimate")

# forplot <- rbind(gradPath, coefPath)
# ggplot(forplot) +
#   geom_line(aes(x = Iteration, y = Gradient, col = factor(Variable))) +
#   theme_bw() + geom_hline(yintercept = 0) +
#   facet_wrap(~ Value, scales = "free") +
#   theme(legend.position = "none") + ylab("")
# ggsave(filen = "vartex/gradientPlot.pdf", width = 7, height = 3)

naivefit <- lm(y ~ X[, select] - 1)
sds <- sqrt(diag(vcov(naivefit)))
naiveCI <- matrix(nrow = sum(select), ncol = 2)
naiveCI[, 1] <- coef(naivefit) - 2 * sds
naiveCI[, 2] <- coef(naivefit) + 2 * sds

trueCoef <- coef(lm(mu ~ X[, select] - 1))
plotEstimates(estimates = list(naive = fit$naive, mle = fit$mle, true = trueCoef),
              CIs = list(naive = naiveCI, bootstrap = fit$naiveBootCI, poly = fit$polyCI))
# ggsave(filen = "vartex/bootCI.pdf", width = 7, height = 3)
ggsave(filen = "vartex/bootCIwPoly.pdf", width = 7, height = 3)

