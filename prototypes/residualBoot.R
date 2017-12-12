# psatThreshold <- function(X, y, ysig, thresholdLevel, test_level, ...) {
#   psatfit <- psatGLM(X, y, test = "wald",
#                       family = "gaussian",
#                       resid_sd = c("null"),
#                       pval_threshold = test_level,
#                       estimate_type = c("naive"),
#                       pvalue_type = c("polyhedral"),
#                       ci_type = c("naive"),
#                       verbose = FALSE)
#   pval <- psatfit$polyPval
#   coef <- psatfit$naive
#   coef[pval < thresholdLevel] <- 0
# }

getCover <- function(ci, truth) {
  cover <- 0
  for(i in 1:length(truth)) {
    if(ci[i, 1] < truth[i] & ci[i, 2] > truth[i]) {
      cover <- cover + 1
    } else {
      cat(i, " ")
    }
  }
  return(cover / length(truth))
}

p <- 20
n <- 50
pthreshold <- 0.01
snr <- 0.1
sparsity <- 2
rho <- 0

X <- matrix(rnorm(n * p), ncol = p)
sigma <- rho^as.matrix(dist(1:p))
sqrtsig <- expm::sqrtm(sigma)
X <- X %*% sqrtsig
X <- scale(X)
true <- rep(0, p)
nonzero <- sample.int(p, sparsity)
true[nonzero] <- (1 - 2 * rbinom(sparsity, 1, 0.5)) * rexp(sparsity)
mu <- as.numeric(X %*% true)
ysig <- sqrt(var(mu) / snr)

XtX <- t(X) %*% X
suffCov <- XtX * ysig^2
suffPrecision <- solve(suffCov)
testPval <- 1
while(testPval > pthreshold) {
  y <- mu + rnorm(n, sd = ysig)
  y <- y - mean(y)
  suffStat <- t(X) %*% y
  waldstat <- as.numeric(t(suffStat) %*% suffPrecision %*% suffStat)
  testPval <- pchisq(waldstat, df = p, lower.tail = FALSE)
  print(testPval)
}

nsamps <- 1000
fit <- residualConditionalBootstrap(X, y, ysig, selected = NULL,
                                    cilevel = 0.05, thresholdLevel = 0.001/p,
                                    bootSamples = nsamps, verbose = TRUE,
                                    precision = XtX / ysig^2, testLevel = pthreshold)
print(nsamps / fit$tries)

naivefit <- lm(y ~ X - 1)
naive <- coef(naivefit)
naiveCov <- solve(t(X) %*% X) * ysig^2
naiveSD <- sqrt(diag(naiveCov))
naiveCI <- matrix(nrow = p, ncol = 2)
naiveCI[, 1] <- naive - naiveSD * qnorm(0.975)
naiveCI[, 2] <- naive + naiveSD * qnorm(0.975)

getCover(fit$ci, true)
getCover(naiveCI, true)

