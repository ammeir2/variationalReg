library(mvtnorm)
library(magrittr)

# Functions ----
checkInclusion <- function(cis, true) {
  inside <- 0
  for(i in 1:nrow(cis)) {
    if(cis[i, 1] < true[i] & cis[i, 2] > true[i]) {
      inside <- inside + 1
    }
  }

  return(inside / nrow(cis))
}

# Generating Data ------
p <- 50
n <- 100
select <- 10
rho <- 0.3
snr <- 0.5
sparsity <- 3
sigma <- rho^as.matrix(dist(1:p))

reps <- 100
results <- matrix(nrow = reps, ncol = 4)
names(results) <- c("poly", "nboot", "mle", "ostep")
set.seed(1)
for(m in 1:reps) {
  true <- rep(0, p)
  true[sample.int(p, sparsity)] <- rnorm(sparsity)
  X <- rmvnorm(n, sigma =  sigma)
  X <- scale(X)
  mu <- as.numeric(X %*% true)
  ysig <- sqrt(var(mu) / snr)
  y <- rnorm(n, mean = mu, sd = ysig)
  y <- y - mean(y)

  suffStat <- as.numeric(t(X) %*% y)
  selected <- (1:p) %in% order(abs(suffStat), decreasing = TRUE)[1:select]
  threshold <- (max(abs(suffStat[!selected])) + min(abs(suffStat[selected]))) / 2
  threshold <- max(abs(suffStat[!selected]))

  # Analysis ------
  Xm <- X[, selected]
  trueProj <- lm(mu ~ Xm - 1) %>% coef()
  fit <- approxConditionalMLE(X, y, ysig, threshold,
                              thresholdLevel = 0.01, cilevel = 0.05,
                              varCI = FALSE,
                              bootSamples = 2000,
                              verbose = TRUE,
                              computeMLE = TRUE)

  # Experimenting with one step approximation -------
  weirdCI <- oneStepCI(fit)
  # trueCoef <- round(trueProj, 4)
  mle <- round(fit$mle, 4)
  # naive <- lm(y ~ X[, selected] - 1) %>% coef()
  # tCoef <- fit$sampCoef
  #
  # condSamp <- rbind(suffStat, fit$suffSamp)
  # suffsamp <- condSamp / nrow(X)
  # forQuantiles <- apply(suffsamp, 2, function(x) x - mean(x))
  # variance <- var(sqrt(nrow(X)) * forQuantiles)
  # A <- diag(p)
  # A[selected, ] <- variance[selected, ]
  # Ainv <- solve(A)
  # forQuantiles <- forQuantiles %*% Ainv * ysig^2
  # obs <- forQuantiles[1, ]
  # forQuantiles <- forQuantiles[-1, ]
  # cilevel <- 0.05
  # quantiles <- apply(forQuantiles[, selected, drop = FALSE], 2,
  #                    function(x) quantile(x, probs = c(1 - cilevel / 2, cilevel / 2)))
  # oneStep <- tCoef + obs[selected]
  # weirdCI <- matrix(nrow = sum(selected), ncol = 2)
  # for(i in 1:nrow(weirdCI)) {
  #   weirdCI[i, ] <- oneStep[i] - quantiles[, i]
  # }

  # cbind(mle, oneStepMLE, naive)
  # cbind(fit$polyCI, trueProj)
  # cbind(fit$naiveBootCI, trueProj)
  # cbind(fit$mleCI, trueProj)
  # cbind(weirdCI, trueProj)
  poly <- checkInclusion(fit$polyCI, trueProj)
  boot <- checkInclusion(fit$naiveBootCI, trueProj)
  mle <- checkInclusion(fit$mleCI, trueProj)
  ostep <- checkInclusion(weirdCI, trueProj)
  results[m, ] <- c(poly, boot, mle, ostep)
  print(colMeans(results[1:m, , drop = FALSE]))
}


