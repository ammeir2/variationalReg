library(variationalReg)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
seed <- as.numeric(seed)

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

    # Estimating --------------
    fit <- NULL
    try(fit <- approxConditionalMLE(X, y, ysig, threshold, thresholdLevel = 0.01 / nselect,
                                verbose = FALSE, bootSamples = 2000,
                                #true = true, trueCoef = projTrue,
                                varCI = TRUE))
    # fit$varBootCI <- fit$naiveBootCI
    polyCI <- fit$polyCI
    mle <- NULL
    try(mle <- exactMSmle(X, y, ysig, threshold, nsteps = 2000, stepCoef = 0.01, stepRate = 0.6,
                      verbose = FALSE))
    if(is.null(fit) | is.null(mle)) {
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
    hybrid <- summary(fit, ci_types = c("poly", "naiveBoot"), ci_level = 0.05 / 2)
    hybrid <- cbind(pmax(hybrid[[1]][, 1], hybrid[[2]][, 1]), pmin(hybrid[[1]][, 2], hybrid[[2]][, 2]))

    # Reporting -----------------
    estimates <- data.frame(naive = fit$naive, mle = mle$mle, var = fit$mEst,
                            true = true[selected], projTrue = projTrue)
    cis <- list(naive = naivenaiveCI,
                naiveBoot = fit$naiveBootCI, varBoot = fit$varBootCI,
                hybrid = hybrid, poly = polyCI)
    results[[m]] <- list(config = config, estimate = estimates, cis = cis)

    # Evaluating ----------------
    print(c(m = m, config))
    coverage <- coverage * (m - 1) / m + sapply(cis, getCover, projTrue) / m
    print(coverage)
    mse <- mse * (m - 1)/m  + apply(estimates[, -(4:5)], 2, function(x) sqrt(mean((x - projTrue)^2)))/m
    print(mse)
    print(sapply(cis, function(x) mean(x[, 2] - x[, 1])))
  }

  return(results)
}

configurations <- expand.grid(n = c(200, 400, 800, 1600),
                              p = c(100),
                              snr = c(0.05, 0.2, 0.8),
                              sparsity = c(1, 2, 4),
                              covtype = c(2),
                              rho = c(0, 0.35, 0.7),
                              nselect = c(10),
                              reps = 1)
set.seed(seed)
subconfig <- configurations[sample.int(nrow(configurations), 20), ]
results <- apply(subconfig, 1, run.sim)
filename <- paste("results/variationalSim_F", seed, ".rds", sep = "")
saveRDS(object = results, file = filename)




