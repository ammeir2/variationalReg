# library(variationalReg)
# args <- commandArgs(TRUE)
# eval(parse(text=args[[1]]))
# seed <- as.numeric(seed)

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

run.sim <- function(config) {
  # Getting variables ---------
  p <- config[["p"]]
  n <- config[["n"]]
  nselect <- config[["nselect"]]
  snr <- config[["snr"]]
  sparsity <- config[["sparsity"]]
  reps <- config[["reps"]]
  covtype <- config[["covtype"]]
  rho <- config[["rho"]]

  # Setting up X covariance --------
  if(covtype == 1) {
    sqrtsig <- diag(p)
  } else if(covtype == 2) {
    sigma <- rho^as.matrix(dist(1:p))
    sqrtsig <- expm::sqrtm(sigma)
  } else if(covtype == 3) {
    rho <- 0.3
    sigma <- matrix(rho, nrow = p, ncol = p)
    diag(sigma) <- 1
    sqrtsig <- expm::sqrtm(sigma)
  }

  results <- vector("list", reps)
  coverage <- 0
  mse <- 0
  m <- 0
  while(m < reps) {
    m <- m + 1
    print(c(m = m, config))
    # Generating Data ------------
    X <- matrix(rnorm(n * p), ncol = p)
    X <- X %*% sqrtsig
    X <- scale(X)
    true <- rep(0, p)
    nonzero <- sample.int(p, sparsity)
    true[nonzero] <- (1 - 2 * rbinom(sparsity, 1, 0.5)) * rexp(sparsity)
    true <- true / sum(abs(true))
    trueCoef <- true
    mu <- as.numeric(X %*% true)
    ysig <- sqrt(var(mu) / snr)
    y <- mu + rnorm(n, sd = ysig)
    y <- y - mean(y)
    suffStat <- t(X) %*% y
    threshold <- mean(sort(abs(suffStat), decreasing = TRUE)[nselect:(nselect + 1)])
    selected <- abs(suffStat) > threshold
    projTrue <- round(coef(lm(mu ~ X[, selected] - 1)), 6)
    ysplit <- (mu + rnorm(n, sd = ysig)) - mean(y)

    # Estimating --------------
    fit <- NULL
    try(fit <- approxConditionalMLE(X, y, ysig, threshold, thresholdLevel = 0.01 / nselect,
                                verbose = TRUE, bootSamples = 2000,
                                thresholdContrast = TRUE,
                                #true = true, trueCoef = projTrue,
                                varCI = TRUE))
    # fit$varBootCI <- fit$naiveBootCI
    polyCI <- fit$polyCI
    polyBoot <- fit$polyBootCI
    mle <- NULL
    if(is.null(fit)) {
      m <- m - 1
      next
    }

    # Naive CI ---------
    naiveFit <- lm(y ~ X[, selected] - 1)
    naiveSD <- sqrt(diag(solve(t(X[, selected]) %*% X[, selected]))) * ysig
    naivenaiveCI <- matrix(nrow = nselect, ncol = 2)
    naivenaiveCI[, 1] <- coef(naiveFit) - qnorm(.975) * naiveSD
    naivenaiveCI[, 2] <- coef(naiveFit) + qnorm(.975) * naiveSD

    # Computing hybrid CI -------
    # hybrid <- summary(fit, ci_types = c("poly", "naiveBoot"), ci_level = 0.05 / 2)
    # hybrid <- cbind(pmax(hybrid[[1]][, 1], hybrid[[2]][, 1]), pmin(hybrid[[1]][, 2], hybrid[[2]][, 2]))

    # Split estimator ------
    indfit <- lm(ysplit ~ X[, selected] - 1)
    indCoef <- coef(indfit)
    indSD <- sqrt(diag(vcov(indfit)))
    indCI <- matrix(nrow = nselect, ncol = 2)
    indCI[, 1] <- indCoef - qnorm(.975) * indSD
    indCI[, 2] <- indCoef + qnorm(.975) * indSD

    # Reporting -----------------
    estimates <- data.frame(naive = fit$naive, mle = fit$mle, var = fit$mEst,
                            ind = indCoef,
                            true = true[selected], projTrue = projTrue)
    cis <- list(naive = naivenaiveCI,
                naiveBoot = fit$naiveBootCI, varBoot = fit$varBootCI,
                poly = fit$polyCI,
                #polyBoot = fit$polyBoot,
                ind = indCI)
    results[[m]] <- list(config = config, estimate = estimates, cis = cis)

    # Evaluating ----------------
    coverage <- coverage * (m - 1) / m + sapply(cis, getCover, projTrue) / m
    print(coverage)
    mse <- mse * (m - 1)/m  + apply(estimates[, -(5:6)], 2, function(x) sqrt(mean((x - projTrue)^2)))/m
    print(mse)
    print(sapply(cis, function(x) mean(x[, 2] - x[, 1])))
  }

  return(results)
}

configurations <- expand.grid(n = 200,
                              p = c(100),
                              snr = 2^((-10):1),
                              # snr = c(0.25, 0.5, 1, 2),
                              sparsity = c(1, 2, 4, 8),
                              covtype = c(2),
                              rho = c(0, 0.35, 0.7),
                              # rho = c(0.7),
                              nselect = c(10),
                              reps = 1)
set.seed(seed)
runif(4)
subconfig <- configurations[sample.int(nrow(configurations), 40), ]
results <- apply(subconfig, 1, run.sim)
filename <- paste("results/variationalSim_univZ_A", seed, ".rds", sep = "")
saveRDS(object = results, file = filename)




