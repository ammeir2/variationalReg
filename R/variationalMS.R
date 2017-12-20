approxConditionalMLE <- function(X, y, ysig, threshold,
                                 thresholdMethod = c("poly", "naiveAdj"),
                                 thresholdLevel = 0.01, cilevel = 0.05,
                                 true = NULL, trueCoef = NULL, varCI = TRUE,
                                 bootSamples = 2000,
                                 verbose = TRUE, thresholdContrast = FALSE) {
  # Variational Estimate --------------------
  p <- ncol(X)
  suffStat <- t(X) %*% y
  selected <- abs(suffStat) > threshold
  Xs <- X[, selected]
  suffCov <- t(X) %*% X * ysig^2
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
  selectedCov <- XtX * ysig^2
  selectedPrecision <- XtXinv / ysig^2
  suffEst <- msVariationalOptim(Xsy, selectedPrecision, selectedCov, threshold)
  mvarEst <- as.numeric(XtXinv %*% suffEst)

  # Computing conditionalMLE
  mle <- exactMSmle(X, y, ysig, threshold,
                    nsteps = 1000, stepCoef = 0.01, stepRate = 0.6,
                    verbose = verbose)$mle

  # Computing Variational CI ------------------
  naivesd <- summary(naiveFit)$coefficients[, 2]
  hardthreshold <- naive - sign(naive) * naivesd * qnorm(1 - thresholdLevel)
  thresholdCoef <- mvarEst
  thresholdCoef[sign(naive) != sign(hardthreshold)] <- 0
  if(is.null(trueCoef)) {
    thresholdMethod <- thresholdMethod[1]
    thresholdMean <- rep(0, p)
    if(thresholdMethod == "poly") {
      polyfit <- polyhedralMS(y, X, as.numeric(suffStat),
                              suffCov, ysig, selected, Eta = NULL,
                              delta = 10^-4,
                              computeCI = TRUE, computeBootCI = FALSE,
                              level = 1 - cilevel)
      if(!thresholdContrast) {
        thresholdCoef <- mle
        polypval <- polyfit$pval
        thresholdCoef[polypval > thresholdLevel] <- 0
      } else {
        mleMean <- as.numeric(XtX %*% mle)
        Eta <- matrix(0, ncol = sum(selected), nrow = p)
        Eta[selected, ] <- diag(sum(selected))
        # contrastpval <- polyhedralMS(y, X, as.numeric(suffStat),
        #                              suffCov, ysig, selected, Eta = Eta,
        #                              computeCI = FALSE, computeBootCI = FALSE)$pval
        contrastpval <- univariateScreenPval(suffStat, suffCov, selected, threshold)
        contrastpval[is.nan(contrastpval)] <- 0
        mleMean[contrastpval > thresholdLevel] <- 0
        # print(mean(contrastpval < thresholdLevel))
        thresholdCoef <- XtXinv %*% mleMean
      }
    } else if(thresholdMethod == "naiveAdj") {
      hardthreshold <- naive - sign(naive) * naivesd * qnorm(1 - thresholdLevel)
      thresholdCoef <- mvarEst
      thresholdCoef[sign(naive) != sign(hardthreshold)] <- 0
    } else {
      stop("thresholdMethod not supported!")
    }
    thresholdMean[selected] <- XtX %*% thresholdCoef
  } else {
    thresholdMean <- as.numeric(suffCov %*% true) / ysig^2
  }
  sampthreshold <- matrix(threshold, nrow = p, ncol = 2)
  sampthreshold[, 1] <- -sampthreshold[, 1]
  diagvar <- diag(suffCov)
  thresholdMean <- thresholdMean / sqrt(diagvar)
  sampthreshold <- sampthreshold / sqrt(diagvar)
  sampCov <- cov2cor(suffCov)
  sampPrecision <- solve(sampCov)
  if(verbose) print("Sampling!")
  restarts <- sum(selected) * 2
  samples <- vector('list', restarts)
  if(verbose) pb <- txtProgressBar(min = 0, max = restarts, style = 3)
  for(i in 1:restarts) {
    start <- rep(0, p)
    samporder <- (1:p)[order(runif(p))]
    condSamp <- mvtSampler(y = start, mu = as.numeric(thresholdMean),
                           selected = as.integer(selected),
                           threshold = sampthreshold,
                           precision = sampPrecision, nsamp = ceiling(bootSamples / restarts),
                           burnin = 200, trim = 40, samporder = samporder, verbose = FALSE)
    samples[[i]] <- condSamp
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)
  condSamp <- do.call("rbind", samples)
  # print(apply(condSamp[, selected], 2, function(x) min(mean(x < 0), mean(x > 0))))
  for(i in 1:ncol(condSamp)) {
    condSamp[, i] <- condSamp[, i] * sqrt(diagvar[i])
  }

  # Naive Cond Boot -----------------
  suffSamp <- condSamp[, selected]
  naiveBoot <- suffSamp %*% XtXinv
  if(is.null(true)) {
    if(thresholdContrast) {
      # thresholdedNaive <- suffStat[selected]
      # thresholdedNaive[mleMean == 0] <- 0
      # thresholdedNaive <- XtXinv %*% thresholdedNaive
      thresholdedNaive <- thresholdCoef
    } else {
      thresholdedNaive <- naive
      thresholdedNaive[thresholdCoef == 0] <- 0
    }
  } else {
    thresholdedNaive <- trueCoef
  }
  naiveBootCI <- computeBootCI(naiveBoot, thresholdedNaive, naive, cilevel)

  # Variational Cond Boot ----------------
  if(varCI) {
    # mvarCoef <- mvarEst
    # mvarCoef[thresholdCoef == 0] <- 0
    varBoot <- computeVarBootSample(suffSamp, selectedPrecision, selectedCov, XtXinv, threshold, verbose)
    varBootCI <- computeBootCI(varBoot, thresholdCoef, mvarEst, cilevel)
  } else {
    varBootCI <- NULL
    varBoot <- NULL
  }

  polyCI <- polyfit$polyCI
  polyBootCI <- polyfit$polyBootCI

  result <- list(mle = mle, mEst = mvarEst, naive = naive,
                 naiveBootCI = naiveBootCI,
                 varBootCI = varBootCI,
                 polyCI = polyfit$ci,
                 polyBootCI = polyfit$bootCI,
                 naiveBoot = naiveBoot,
                 varBoot = varBoot,
                 sampCoef = thresholdCoef,
                 cilevel = cilevel,
                 ysig = ysig,
                 suffStat = suffStat,
                 suffCov = suffCov,
                 threshold = threshold,
                 selected = selected,
                 call = match.call())
  class(result) <- "varMS"
  return(result)
}

computeBootCI <- function(sample, sampparam, estimate, level) {
  ci <- matrix(nrow = ncol(sample), ncol = 2)
  for(i in 1:nrow(ci)) {
    ci[i, ] <- estimate[i] - quantile(sample[, i] - sampparam[i], c(1 - level / 2, level / 2))
  }
  return(ci)
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
  # dens <- sum(dnorm(y, yhat, ysig, log = TRUE))
  suffMu <- as.numeric(suffCov %*% param) / ysig^2
  sds <- sqrt(diag(suffCov))
  dens <- sum(dnorm(as.numeric(t(Xs) %*% y), suffMu, sds, log = TRUE))
  prob <- pnorm(-threshold, suffMu, sds) + pnorm(threshold, suffMu, sds, lower.tail = FALSE)
  return(-dens + sum(log(prob)))
}

computeVarBootSample <- function(suffSamp, suffPrecision, suffCov, XtXinv, threshold, verbose) {
  if(verbose) {
    print("Boostrapping variational estimate!")
    pb <- txtProgressBar(min = 0, max = nrow(suffSamp), style = 3)
  }
  varBoot <- matrix(nrow = nrow(suffSamp), ncol = ncol(XtXinv))
  for(i in 1:nrow(suffSamp)) {
    if(verbose) setTxtProgressBar(pb, i)
    varBoot[i, ] <- msVariationalOptim(suffSamp[i, ], suffPrecision, suffCov, threshold)
  }
  if(verbose) close(pb)
  varBoot <- varBoot %*% XtXinv
  return(varBoot)
}

msVariationalOptim <- function(Xsy, suffPrecision, suffCov, threshold) {
  optim(par = Xsy, fn = compSuffMeanLoglikM,
        method = "Nelder-Mead",
        suffStat = Xsy,
        suffPrecision = suffPrecision,
        suffCov = suffCov,
        threshold = threshold)$par
}

computeNaiveCI <- function(suffStat, suffCov, selected, ysig, level) {
  XtX <- suffCov[selected == 1, ][, selected == 1] / ysig^2
  XtXinv <- solve(XtX)
  coef <- as.numeric(XtXinv %*% suffStat[selected])
  coefcov <- XtXinv * ysig^2
  quant <- qnorm(1 - level, sd = sqrt(diag(coefcov)))
  ci <- matrix(nrow = nrow(XtX), ncol = 2)
  ci[, 1] <- coef - quant * diag(coefcov)
  ci[, 2] <- coef + quant * diag(coefcov)
  return(ci)
}


univariateScreenPval <- function(suffStat, suffCov, selected, threshold) {
  y <- suffStat[selected]
  sds <- sqrt(diag(suffCov))[selected]
  pselect <- pnorm(-threshold, mean = 0, sd = sds, lower.tail = TRUE)
  pselect <- pselect + pnorm(threshold, mean = 0, sd = sds, lower.tail = FALSE)
  prange <- 2* pnorm(-abs(y), mean = 0, sd = sds)
  return(prange / pselect)
}





