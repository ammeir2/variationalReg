approxConditionalMLE <- function(X, y, ysig, threshold, thresholdLevel = 0.01,
                                 true = NULL, trueCoef = NULL, varCI = TRUE,
                                 bootSamples = 1000,
                                 verbose = TRUE) {
  # Variational Estimate --------------------
  p <- ncol(X)
  suffStat <- t(X) %*% y
  selected <- abs(suffStat) > threshold
  Xs <- X[, selected]
  suffCov <- t(Xs) %*% Xs * ysig^2
  XtX <- t(Xs) %*% Xs
  XtXinv <- solve(XtX)

  naiveFit <- lm(y ~ Xs - 1)
  naive <- coef(naiveFit)
  suffPrecision <- solve(suffCov)

  # mvarEst <- optim(par = naive, fn = compMloglik,
  #                  method = "Nelder-Mead",
  #                  y = y, Xs = Xs, ysig = ysig, suffCov = suffCov,
  #                  threshold = threshold)

  Xsy <- as.numeric(suffStat[selected, ])
  suffEst <- msVariationalOptim(Xsy, suffPrecision, suffCov, threshold)
  mvarEst <- as.numeric(XtXinv %*% suffEst)

  # Computing Variational CI ------------------
  naivesd <- summary(naiveFit)$coefficients[, 2]
  hardthreshold <- naive - sign(naive) * naivesd * qnorm(1 - thresholdLevel)
  thresholdCoef <- naive
  thresholdCoef[sign(naive) != sign(hardthreshold)] <- 0
  if(is.null(trueCoef)) {
    thresholdMean <- rep(0, p)
    thresholdMean[selected] <- XtX %*% thresholdCoef
  } else {
    thresholdMean <- as.numeric(t(X) %*% X %*% true)
  }
  sampthreshold <- matrix(threshold, nrow = p, ncol = 2)
  sampthreshold[, 1] <- -sampthreshold[, 1]
  sampCov <- t(X) %*% X * ysig^2
  diagvar <- diag(sampCov)
  thresholdMean <- thresholdMean / sqrt(diagvar)
  sampthreshold <- sampthreshold / sqrt(diagvar)
  sampCov <- cov2cor(sampCov)
  sampPrecision <- solve(sampCov)
  if(verbose) print("Sampling!")
  condSamp <- mvtSampler(y = as.numeric(suffStat / sqrt(diagvar)), mu = as.numeric(thresholdMean),
                         selected = as.integer(selected),
                         threshold = sampthreshold,
                         precision = sampPrecision, nsamp = bootSamples,
                         burnin = 1000, trim = 10, verbose = verbose)
  for(i in 1:ncol(condSamp)) {
    condSamp[, i] <- condSamp[, i] * sqrt(diagvar[i])
  }

  # Naive Cond Boot -----------------
  suffSamp <- condSamp[, selected]
  naiveBoot <- suffSamp %*% XtXinv
  naiveBootCI <- matrix(nrow = ncol(Xs), ncol = 2)
  if(is.null(true)) {
    thresholdedNaive <- naive
    thresholdedNaive[thresholdCoef == 0] <- 0
  } else {
    thresholdedNaive <- trueCoef
  }
  for(i in 1:nrow(naiveBootCI)) {
    naiveBootCI[i, ] <- naive[i] - quantile(naiveBoot[, i] - thresholdedNaive[i], c(0.975, 0.025))
  }

  # Variational Cond Boot ----------------
  if(varCI) {
    if(verbose) {
      print("Boostrapping variational estimate!")
      pb <- txtProgressBar(min = 0, max = nrow(suffSamp), style = 3)
    }
    varBoot <- matrix(nrow = nrow(suffSamp), ncol = sum(selected))
    for(i in 1:nrow(suffSamp)) {
      if(verbose) setTxtProgressBar(pb, i)
      varBoot[i, ] <- msVariationalOptim(suffSamp[i, ], suffPrecision, suffCov, threshold)
    }
    varBoot <- varBoot %*% XtXinv
    varBootCI <- matrix(nrow = ncol(Xs), ncol = 2)
    if(!is.null(true)) {
      thresholdCoef <- trueCoef
    }

    for(i in 1:nrow(varBootCI)) {
      varBootCI[i, ] <- mvarEst[i] - quantile(varBoot[, i] - thresholdCoef[i], c(0.975, 0.025))
    }
  } else {
    varBootCI <- NULL
    varBoot <- NULL
  }

  return(list(mEst = mvarEst, naive = naive,
              naiveBootCI = naiveBootCI,
              varBootCI = varBootCI,
              naiveBoot = naiveBoot,
              varBoot = varBoot))
}

compSuffMeanLoglikM <- function(mean, suffStat, suffPrecision, suffCov, threshold) {
  diff <- suffStat - mean
  dens <- -0.5 * t(diff) %*% suffPrecision %*% diff
  sds <- sqrt(diag(suffCov))
  prob <- pnorm(-threshold, mean, sds) + pnorm(threshold, mean, sds, lower.tail = FALSE)
  return(-dens + sum(log(prob)))
}

compMloglik <- function(param, y, Xs, ysig, suffCov, threshold) {
  yhat <- as.numeric(Xs %*% param)
  dens <- sum(dnorm(y, yhat, ysig, log = TRUE))

  suffMu <- as.numeric(suffCov %*% param) / ysig^2
  sds <- sqrt(diag(suffCov))
  prob <- pnorm(-threshold, suffMu, sds) + pnorm(threshold, suffMu, sds, lower.tail = FALSE)
  return(-dens + sum(log(prob)))
}

msVariationalOptim <- function(Xsy, suffPrecision, suffCov, threshold) {
  optim(par = Xsy, fn = compSuffMeanLoglikM,
        method = "Nelder-Mead",
        suffStat = Xsy,
        suffPrecision = suffPrecision,
        suffCov = suffCov,
        threshold = threshold)$par
}
