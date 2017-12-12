# args <- commandArgs(TRUE)
# eval(parse(text=args[[1]]))
# seed <- as.numeric(seed)

leapsVerify <- function(X, y, ysig, selected, ...) {
  lfit <- leaps(x = X, y = y, int = FALSE)
  best <- lfit$which[which.min(lfit$Cp), ]
  # print(as.numeric(which(best)))
  # print(as.numeric((which(selected))))
  # print(all(selected == best))
  return(all(selected == best))
}

library(PSAT)
library(variationalReg)
library(leaps)

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
  pthreshold <- config[["pthreshold"]]
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
  sizes <- 0
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
    mu <- as.numeric(X %*% true)
    ysig <- sqrt(var(mu) / snr)

    X <- matrix(rnorm(n * p), ncol = p)
    X <- X %*% sqrtsig
    X <- scale(X)
    true <- rep(0, p)
    nonzero <- sample.int(p, sparsity)
    true[nonzero] <- (1 - 2 * rbinom(sparsity, 1, 0.5)) * rexp(sparsity)
    mu <- as.numeric(X %*% true)
    ysig <- sqrt(var(mu) / snr)
    y <- mu + rnorm(n, sd = ysig)
    y <- y - mean(y)

    XtX <- t(X) %*% X
    suffCov <- XtX * ysig^2
    suffPrecision <- solve(suffCov)

    lfit <- leaps(x = X, y = y, int = FALSE)
    selected <- lfit$which[which.min(lfit$Cp), ]

    # Estimating --------------
    fit <- NULL
    nsamps <- 500
    precision <- solve(t(X) %*% X * ysig^2)
    which(selected)
    fit <- residualConditionalBootstrap(X, y, ysig, selected = selected,
                                        cilevel = 0.05, thresholdLevel = pthreshold^2/p,
                                        bootSamples = nsamps, verbose = TRUE,
                                        selectionF = leapsVerify,
                                        testLevel = pthreshold)
    print(nsamps / fit$tries)

    # Naive -----
    naivefit <- lm(y ~ X - 1)
    naive <- coef(naivefit)
    naiveCov <- solve(t(X) %*% X) * ysig^2
    naiveSD <- sqrt(diag(naiveCov))
    naiveCI <- matrix(nrow = p, ncol = 2)
    naiveCI[, 1] <- naive - naiveSD * qnorm(0.975)
    naiveCI[, 2] <- naive + naiveSD * qnorm(0.975)

    # PSAT ---------------
    psat <- PSAT::psatGLM(X, y, test = "wald",
                          resid_sd = ysig,
                          pval_threshold = pthreshold,
                          verbose = FALSE)

    # Reporting ---------
    cis <- list(naive = naiveCI,
                switch = psat$switchCI,
                poly = psat$polyCI,
                boot = fit$ci)
    results[[m]] <- list(config = config, cis = cis)

    # Evaluating ----------------
    print(c(m = m, config))
    coverage <- coverage * (m - 1) / m + sapply(cis, getCover, true) / m
    print(coverage)
    sizes <- sizes * (m - 1) / m + (sapply(cis, function(x) mean(x[, 2] - x[, 1]))) / m
    print(sizes)
  }

  return(results)
}

configurations <- expand.grid(n = c(100),
                              p = c(8),
                              snr = c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8),
                              sparsity = c(1, 2, 4),
                              covtype = c(2),
                              rho = c(0.35),
                              pthreshold = c(0.1, 0.01, 0.001),
                              reps = 1)
set.seed(seed)
results <- apply(configurations, 1, run.sim)
filename <- paste("results/variationalSimAgg_A", seed, ".rds", sep = "")
saveRDS(object = results, file = filename)



