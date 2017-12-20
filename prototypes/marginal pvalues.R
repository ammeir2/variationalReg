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

univpval <- function(y, mu, sd, threshold, cilevel = NULL) {
  denom <- pnorm(-threshold, mean = mu, sd = sd) + pnorm(threshold, mean = mu, sd = sd, lower.tail = FALSE)
  numerator <- pnorm(min(y, -threshold), mean = mu, sd = sd)
  if(y > threshold) {
    numerator <- numerator - pnorm(threshold, mean = mu, sd = sd) + pnorm(y, mean = mu, sd = sd)
  }
  pval <- numerator / denom
  pval <- 2 * min(pval, 1 - pval)
  if(is.null(cilevel)) {
    return(pval)
  } else {
    return(pval > 1 - cilevel)
  }
}

run.sim <- function(config) {
  # Getting variables ---------
  p <- config[["p"]]
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
  error <- 0
  pow <- 0
  cover <- 0
  m <- 0
  mt <- 0
  mp <- 0
  while(m < reps) {
    m <- m + 1
    # Generating Data ------------
    y <- rnorm(p)
    y <- as.numeric(sqrtsig %*% y)
    nonzero <- sample.int(p, sparsity)
    if(sparsity > 1) {
      signal <-  (1 - 2 * rbinom(sparsity, 1, 0.5)) * rexp(sparsity)
      signal <- signal - mean(signal)
      signal <- signal / sd(signal) * sqrt(snr)
    } else {
      signal <- snr
    }
    mu <- rep(0, p)
    mu[nonzero] <- signal
    y <- y + mu
    threshold <- qnorm(1 - 0.05 / 2, sd = 1)
    selected <- abs(y) > threshold
    if(all(!selected)){
      m <- m - 1
      next
    }

    # Marginal pvals --------------
    true <- mu[selected]
    whichSelected <- which(selected)
    adjpvals <- numeric(sum(selected))
    adjcover <- numeric(sum(selected))
    for(i in 1:sum(selected)) {
      ind <- whichSelected[i]
      adjpvals[i] <- univpval(y[ind], mu = 0, sd = sqrt(sigma[ind, ind]), threshold = threshold)
      adjcover[i] <- univpval(y[ind], mu = true[i], sd = sqrt(sigma[ind, ind]), threshold = threshold, cilevel = 0.95)
    }

    # Naive CI ---------
    naivepvals <- 2 * pnorm(-abs(y[selected]), sd = sqrt(diag(sigma))[selected])
    naivecover <- 2 * pnorm(-abs(y[selected] - true), sd = sqrt(diag(sigma))[selected]) > 0.05

    # Reporting -----------------
    results[[m]] <- list(true = mu[selected], adjpvals = adjpvals, naivepvals = naivepvals)

    # Evaluating ----------------
    bonf <- 0.05
    if(any(true == 0)) {
      mt <- mt + 1
      typeI <- c(adj = mean(adjpvals[true == 0] < bonf),
                 naive = mean(naivepvals[true == 0] < bonf))
      error <- error * (mt - 1) / mt + typeI / mt
    }

    if(any(true != 0)) {
      mp <- mp + 1
      power <- c(adj = mean(adjpvals[true != 0] < bonf),
                 naive = mean(naivepvals[true != 0] < bonf))
      pow <- pow * (mp - 1) / mp  + power / mp
    }

    itercover <- c(mean(adjcover), mean(naivecover))
    cover <- cover * (m - 1) / m + itercover / m
  }

  print(c(m = m, config))
  cat("typeI:", error, "\n")
  cat("power:", pow, "\n")
  cat("cover:", cover, "\n")


  return(results)
}

configurations <- expand.grid(p = c(30),
                              snr = c(8, 4, 2, 1, 0.5, 0),
                              sparsity = c(1, 2, 4, 8, 16),
                              covtype = c(2),
                              rho = c(0, 0.4, 0.8),
                              reps = 3 * 10^3)
set.seed(seed)
runif(4)
subconfig <- configurations
results <- apply(subconfig, 1, run.sim)

