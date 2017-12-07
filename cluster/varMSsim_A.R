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

  # Setting up X covariance --------
  if(covtype == 1) {
    sqrtsig <- diag(p)
  } else if(covtype == 2) {
    rho <- 0.8
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
  for(m in 1:reps) {
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
    fit <- approxConditionalMLE(X, y, ysig, threshold, thresholdLevel = 0.01 / nselect,
                                verbose = FALSE)
    polyCI <- polyhedralMS(X, y, ysig, selected, Eta = NULL)
    mle <- exactMSmle(X, y, ysig, threshold, nsteps = 4000, stepCoef = 0.01, stepRate = 0.6,
                      verbose = FALSE)

    naiveFit <- lm(y ~ X[, selected] - 1)
    naivenaiveCI <- matrix(nrow = nselect, ncol = 2)
    naivenaiveCI[, 1] <- coef(naiveFit) - qnorm(.975) * sqrt(diag(vcov(naiveFit)))
    naivenaiveCI[, 2] <- coef(naiveFit) + qnorm(.975) * sqrt(diag(vcov(naiveFit)))

    # Reporting -----------------
    estimates <- data.frame(naive = fit$naive, mle = mle$mle, var = fit$mEst,
                            true = true[selected], projTrue = projTrue)
    cis <- list(naive = naivenaiveCI,
                naiveBoot = fit$naiveBootCI, varBoot = fit$varBootCI,
                poly = polyCI)
    results[[m]] <- list(config = config, estimate = estimates, cis = cis)

    # Evaluating ----------------
    print(c(m = m, config))
    coverage <- coverage * (m - 1) / m + sapply(cis, getCover, projTrue) / m
    print(coverage)
    mse <- mse * (m - 1)/m  + apply(estimates[, -(4:5)], 2, function(x) sqrt(mean((x - projTrue)^2)))/m
    print(mse)
  }

  return(results)
}

configurations <- expand.grid(n = c(250, 500, 1000),
                              p = c(200),
                              snr = c(0.05, 0.2, 0.8),
                              sparsity = c(1, 2, 4),
                              covtype = 1,
                              nselect = c(20),
                              reps = 2)

set.seed(seed)
results <- apply(configurations, 1, run.sim)
filename <- paste("results/variationalSim_B", seed, ".rds", sep = "")
saveRDS(object = results, file = filename)




