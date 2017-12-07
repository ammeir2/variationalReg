exactMSmle <- function(X, y, ysig, threshold,
                       nsteps = 2000, nsamps = 4, stepCoef = 0.02,
                       stepRate = 0.6,
                       meanMethod = c("plugin", "zero"),
                       verbose = TRUE) {
  meanMethod <- meanMethod[1]

  p <- ncol(X)
  suffStat <- t(X) %*% y
  selected <- abs(suffStat) > threshold
  Xs <- X[, selected]
  suffCov <- t(X) %*% X * ysig^2
  XtX <- t(Xs) %*% Xs
  XtXinv <- solve(XtX)

  naiveFit <- lm(y ~ Xs - 1)
  naive <- coef(naiveFit)

  prevGrad <- rep(0, sum(selected))
  betahat <- naive
  if(meanMethod == "plugin") {
    mu <- as.numeric(suffStat)
  } else if(meanMethod == "zero") {
    mu <- rep(0, p)
  } else {
    stop("mean method not supported")
  }

  sampthreshold <- matrix(threshold, nrow = p, ncol = 2)
  sampthreshold[, 1] <- -sampthreshold[, 1]
  prevSamp <- as.numeric(suffStat)
  precision <- solve(suffCov)
  b1 <- 0.9
  b2 <- 0.99
  mt <- 0
  vt <- 0
  betahat <- naive
  estimates <- matrix(nrow = nsteps + 1, ncol = length(betahat))
  estimates[1, ] <- betahat
  if(verbose) {
    print("Computing conditional MLE!")
    pb <- txtProgressBar(min = 0, max = nsteps, style = 3)
  }
  for(m in 1:nsteps) {
    if(verbose) setTxtProgressBar(pb, m)
    mu[selected] <- as.numeric(XtX %*% betahat)
    condSamp <- mvtSampler(y = as.numeric(prevSamp), mu = mu,
                           selected = as.integer(selected),
                           threshold = sampthreshold,
                           precision = precision, nsamp = max(nsamps, 4),
                           burnin = 2, trim = 1, verbose = FALSE)
    prevSamp <- condSamp[nsamps, ]
    grad <- (suffStat[selected] - colMeans(condSamp)[selected])
    mt <- b1 * mt + (1 - b1) * grad
    vt <- b2 * vt + (1 - b2) * grad^2
    betahat <- betahat + (mt / (1 - b1^m)) / (sqrt(vt/(1 - b2^m)) + 10^-6) * stepCoef / m^stepRate

    # correcting signs
    # signs <- sign(naive)
    # betahat <- betahat * signs
    # betahat <- pmax(0, pmin(abs(naive), betahat))
    # betahat <- betahat * signs

    estimates[m + 1, ] <- betahat
  }
  if(verbose) close(pb)

  betahat <- colMeans(estimates[round(nrow(estimates)/2):nrow(estimates), ])
  #print(cbind(betahat, naive))
  return(list(mle = betahat, estimates = estimates, naive = naive))
}
