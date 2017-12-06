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

set.seed(100)
set.seed(102)

p <- 100
n <- 200
pthreshold <- 0.05
nselect <- 20
snr <- 0.2
sparsity <- 3

X <- matrix(rnorm(n * p), ncol = p)
X <- scale(X)
true <- rep(0, p)
nonzero <- sample.int(p, sparsity)
true[nonzero] <- (1 - 2 * rbinom(sparsity, 1, 0.5)) * rexp(sparsity)
trueCoef <- true
mu <- as.numeric(X %*% true)
ysig <- sqrt(var(mu) / snr)
y <- mu + rnorm(n, sd = ysig)
y <- y - mean(y)
suffStat <- t(X) %*% y
threshold <- mean(sort(abs(suffStat), decreasing = TRUE)[nselect:(nselect + 1)])
selected <- abs(suffStat) > threshold
true <- round(coef(lm(mu ~ X[, selected] - 1)), 6)

fit <- approxConditionalMLE(X, y, ysig, threshold, thresholdLevel = 0.1 / nselect /n,
                            #true = trueCoef, trueCoef = true,
                            bootSamples = 400,
                            varCI = TRUE)
polyCI <- polyhedralMS(X, y, ysig, selected, Eta = NULL)
mle <- exactMSmle(X, y, ysig, threshold, nsteps = 2000, stepCoef = 0.05, stepRate = 0.6)
cbind(fit$naive, fit$mEst, true)
naive <- fit$naive
cond <- fit$mEst
adjQuantile <- qnorm(1 - 0.1 / nselect / n)
naiveFit <- lm(y ~ X[, selected] - 1)
adjNaiveCI <- matrix(nrow = nselect, ncol = 2)
adjNaiveCI[, 1] <- coef(naiveFit) - adjQuantile * sqrt(diag(vcov(naiveFit)))
adjNaiveCI[, 2] <- coef(naiveFit) + adjQuantile * sqrt(diag(vcov(naiveFit)))
cbind(adjNaiveCI, true, cond, naive)
cbind(fit$naiveBootCI, true)
cbind(fit$varBootCI, true)
cbind(polyCI, true)
# cbind(adjNaiveCI, true[selected])

naivenaiveCI <- matrix(nrow = nselect, ncol = 2)
naivenaiveCI[, 1] <- coef(naiveFit) - qnorm(.975) * sqrt(diag(vcov(naiveFit)))
naivenaiveCI[, 2] <- coef(naiveFit) + qnorm(.975) * sqrt(diag(vcov(naiveFit)))
getCover(fit$naiveBootCI, true)
getCover(fit$varBootCI, true)
getCover(polyCI, true)
getCover(adjNaiveCI, true)
getCover(naivenaiveCI, true)

mean(polyCI[, 2] - polyCI[, 1])
mean(adjNaiveCI[, 2] - adjNaiveCI[, 1])
mean(fit$naiveBootCI[, 2] - fit$naiveBootCI[, 1])
mean(fit$varBootCI[, 2] - fit$varBootCI[, 1])
mean(2 * sqrt(diag(vcov(naiveFit))) * qnorm(0.975))

sqrt(mean((fit$naive - true)^2))
sqrt(mean((mle$mle - true)^2))
sqrt(mean((fit$mEst - true)^2))
round(cbind(naive = fit$naive, var = fit$mEst, mle = mle$mle, true), 3)

library(ggplot2)
library(reshape2)
library(dplyr)
solutionPath <- mle$estimates
solutionPath <- melt(solutionPath)
names(solutionPath) <- c("iter", "variable", "estimate")
ggplot(solutionPath) + geom_line(aes(x = iter, y = estimate, col = factor(variable), linetype = factor(variable))) +
  geom_hline(yintercept = 0) + theme_bw()

# Plotting estimates/CIs
estimates <- list(naive = naive, variational = fit$mEst, mle = mle$mle, true = true)
CIs <- list(naive = naivenaiveCI, refitBoot = fit$naiveBootCI, varBoot = fit$varBootCI,
            poly = polyCI)
plotEstimates(estimates, CIs) + ylim(min(adjNaiveCI), max(adjNaiveCI))

# Checking RMSE --------------
reps <- 100
naiveRmse <- numeric(reps)
condRmse <- numeric(reps)
condPred <- numeric(reps)
naivePred <- numeric(reps)
for(i in 1:reps) {
  X <- matrix(rnorm(n * p), ncol = p)
  X <- scale(X)
  true <- rep(0, p)
  nonzero <- sample.int(p, sparsity)
  true[nonzero] <- (1 - 2 * rbinom(sparsity, 1, 0.5)) * rexp(sparsity)
  mu <- as.numeric(X %*% true)
  ysig <- sqrt(var(mu) / snr)
  y <- mu + rnorm(n, sd = ysig)
  y <- y - mean(y)
  suffStat <- t(X) %*% y
  threshold <- mean(sort(abs(suffStat), decreasing = TRUE)[nselect:(nselect + 1)])
  selected <- abs(suffStat) > threshold

  fit <- approxConditionalMLE(X, y, ysig, threshold)
  cbind(fit$naive, fit$mEst, true[selected])
  naiveRmse[i] <- sqrt(mean((fit$naive - true[selected])^2))
  condRmse[i] <- sqrt(mean((fit$mEst - true[selected])^2))
  naive <- fit$naive
  cond <- fit$mEst
  condPred[i]<- sqrt(mean((X[, selected] %*% cond - mu)^2))
  naivePred[i] <- sqrt(mean((X[, selected] %*% naive - mu)^2))
  print(c(i,
          mean(log(naiveRmse[1:i]) - log(condRmse[1:i])),
          mean(log(naivePred[1:i]) - log(condPred[1:i]))))
}
cbind(true[selected], naive, fit$mEst)

