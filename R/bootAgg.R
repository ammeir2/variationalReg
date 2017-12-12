waldSelection <- function(X, y, ysig, selected, ...) {
  call <- match.call()
  if(is.null(call$precision)) {
    precision <- solve(t(X) %*% X * ysig^2)
  } else {
    precision <- eval(match.call()$precision)
  }

  if(is.null(call$testLevel)) {
    testLevel <- 0.05
  } else {
    testLevel <- eval(match.call()$testLevel)
  }

  suffStat <- t(X) %*% y
  testStat <- t(suffStat) %*% precision %*% suffStat
  pval <- pchisq(testStat, df = length(suffStat), lower.tail = FALSE)
  #print(c(pval, testLevel))
  return(pval < testLevel)
}

naiveAdjustThreshold <- function(X, y, ysig, selected, thresholdLevel, ...) {
  naivefit <- lm(y ~ X[, selected] - 1)
  coefcov <- vcov(naivefit)
  coef <- coef(naivefit)
  hard <- coef - sign(coef) * sqrt(diag(coefcov)) * qnorm(1 - thresholdLevel)
  coef[sign(hard) != sign(coef)] <- 0
  return(coef)
}

residualConditionalBootstrap <- function(X, y, ysig, selected = NULL,
                                         cilevel = 0.05,
                                         sampCoef = NULL,
                                         thresholdLevel = cilevel^2 / length(y),
                                         selectionF = waldSelection,
                                         thresholdF = naiveAdjustThreshold,
                                         bootSamples = 1000, verbose = TRUE,
                                         ...) {
  if(!selectionF(X, y, ysig, selected, ...)) {
    return(stop("Dataset does not satisfy selection criteria!"))
  }

  if(is.null(selected)) {
    selected <- rep(TRUE, ncol(X))
  }

  Xs <- X[, selected]
  naivefit <- lm(y ~ X[, selected] - 1)
  naive <- coef(naivefit)
  if(is.null(sampCoef)) {
    sampCoef <- thresholdF(X, y, ysig, selected, thresholdLevel, ...)
  } else {
    sampCoef <- sampCoef
  }
  yhat <- as.numeric(Xs %*% sampCoef)
  residuals <- y - yhat
  residuals <- residuals - mean(residuals)

  bootSamp <- matrix(nrow = bootSamples, ncol = length(naive))
  row <- 1
  if(verbose) {
    print("Residual bootstrap")
    pb <- txtProgressBar(min = 0, max = bootSamples, style = 3)
  }

  tries <- 0
  while(row <= bootSamples) {
    tries <- tries + 1
    newy <- yhat + sample(residuals, length(yhat), replace = TRUE)
    if(selectionF(X, newy, ysig, selected, ...)) {
      bootSamp[row, ] <- t(Xs) %*% newy
      row <- row + 1
      if(verbose) setTxtProgressBar(pb, row)
    }
  }
  if(verbose) close(pb)

  XtX <- t(Xs) %*% Xs
  XtXinv <- solve(XtX)
  bootSamp <- bootSamp %*% XtXinv
  bootCI <- computeBootCI(bootSamp, sampCoef, naive, cilevel)
  return(list(ci = bootCI, sample = bootSamp, tries = tries))
}


