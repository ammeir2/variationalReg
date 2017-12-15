polyhedralMS <- function(y, X, suffStat, suffCov, ysig, selected, Eta = NULL, level = 0.95,
                         computeCI = TRUE,
                         computeBootCI = TRUE,
                         verbose = TRUE,
                         delta = 10^-6) {
  p <- length(suffStat)
  sigma <- suffCov

  if(is.null(Eta)) {
    Eta <- matrix(0, ncol = sum(selected), nrow = p)
    Eta[selected, ] <- solve(suffCov[selected, selected] / ysig^2)
  }

  if(computeBootCI) {
    bootEta <- solve(t(X[, selected]) %*% X[, selected]) %*% t(X[, selected])
  }

  # Computing constraints -----------------------
  nconstraints <- 2 * sum(selected) * (p - sum(selected))
  alphaMat <- matrix(nrow = nconstraints, ncol = ncol(Eta))
  Ay <- numeric(nconstraints)
  etavars <- apply(Eta, 2, function(x) t(x) %*% sigma %*% x)
  whichSelected <- which(selected)
  notSelected <- which(!selected)

  constraintInd <- 0
  Arow <- rep(0, ncol(suffCov))
  for(i in 1:sum(selected)) {
    Arow[whichSelected[i]] <- -sign(suffStat[whichSelected][i])
    for(j in 1:sum(!selected)){
      for(s in c(-1, 1)) {
        constraintInd <- constraintInd + 1
        Arow[notSelected[j]] <- s
        Ay[constraintInd] <- t(Arow) %*% suffStat
        alphaMat[constraintInd, ] <- t(Arow) %*% sigma %*% Eta / etavars
      }
      Arow[notSelected[j]] <- 0
    }
    Arow[whichSelected[i]] <- 0
  }

  if(any(Ay > 0)) {
    stop("Selection constraints not satisfied!")
  }

  # Computing confidence intervals -------------------
  if(computeCI) {
    polyCI <- matrix(nrow = ncol(Eta), ncol = 2)
  } else {
    polyCI <- NULL
  }

  if(computeBootCI) {
    bootCI <- matrix(nrow = ncol(Eta), ncol = 2)
    if(verbose) {
      print("Computing polyhedral bootstrap!")
      pb <- txtProgressBar(min = 0, max = ncol(Eta), style = 3)
    }
  } else {
    bootCI <- NULL
  }

  pval <- numeric(ncol(Eta))
  for(i in 1:ncol(Eta)) {
    eta <- Eta[, i]
    etaSigma <- as.numeric(t(eta) %*% sigma %*% eta)
    alpha <- alphaMat[, i]
    theta <- as.numeric(t(eta) %*% suffStat)
    minusSusbset <- alpha < 0
    plusSubset <- alpha > 0
    Vminus <- max((-Ay[minusSusbset] + alpha[minusSusbset] * theta) / alpha[minusSusbset])
    Vplus <- min((-Ay[plusSubset] + alpha[plusSubset] * theta) / alpha[plusSubset])
    if(computeCI) polyCI[i, ] <- findPolyCIlimits(y, theta, eta, etaSigma, Vminus, Vplus, 1 - level, boot = FALSE)
    if(computeBootCI) {
      if(verbose) setTxtProgressBar(pb, i)
      booteta <- as.numeric(bootEta[i, ])
      bootCI[i, ] <- findPolyCIlimits(y, theta, booteta, etaSigma, Vminus, Vplus, 1 - level,
                                      boot = TRUE, delta = delta,
                                      bootsamps = 2000)
    }
    pval[i] <- ptruncnorm(theta, Vminus, Vplus, 0, sqrt(etaSigma))
    pval[i] <- 2 * min(1 - pval[i], pval[i])
  }

  if(computeBootCI & verbose) {
    close(pb)
  }

  return(list(pval = pval, ci = polyCI, bootCI = bootCI))
}


findPolyCIlimits <- function(y, theta, eta, etaSigma, lower, upper, alpha,
                             boot = FALSE, delta = 0, c = 1.0001,
                             bootsamps = 2000) {
  if(boot) {
    center <- mean(y)
    theta <- as.numeric(t(eta) %*% y)
    y <- y - center
    n <- length(y)
    bootsamps <- replicate(bootsamps, sum(eta * sample(y, replace = TRUE)))
  }

  llim <- theta
  steps <- 0
  if(boot) {
    ltempPval <- tibBootQuantile(llim, theta, bootsamps, lower, upper, c, delta = delta)
    step <- sd(bootsamps)
  } else {
    ltempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
    step <- sqrt(etaSigma)
  }
  while(ltempPval < 1 - alpha / 2) {
    steps <- steps + 1
    llim <- llim - step * 0.1
    if(boot) {
      ltempPval <- tibBootQuantile(llim, theta, bootsamps, lower, upper, c, delta = delta)
    } else {
      ltempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
    }
    if(is.nan(ltempPval)) {
      llim <- llim + step * 0.1
      break
    }
    if(steps > 10^5) {
      ltempPval <- NaN
      break
    }
  }

  ulim <- theta
  if(boot) {
    utempPval <- tibBootQuantile(ulim, theta, bootsamps, lower, upper, c, delta = delta)
  } else {
    utempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
  }

  steps <- 0
  while(utempPval > alpha / 2) {
    steps <- steps + 1
    ulim <- ulim + step * 0.1
    if(boot) {
      utempPval <- tibBootQuantile(ulim, theta, bootsamps, lower, upper, c, delta = delta)
    } else {
      utempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
    }
    if(is.nan(utempPval)) {
      ulim <- ulim - step * 0.1
      break
    }
    if(steps > 10^5) {
      utempPval <- NaN
      break
    }
  }

  if(is.nan(utempPval)) {
    uci <- ulim
  } else {
    uci <- NULL
    if(boot) {
      capture.output(invisible(try(uci <- uniroot(f = function(x) tibBootQuantile(x, theta, bootsamps, lower, upper, c, delta = delta) - alpha / 2,
                                                  interval = c(llim, ulim))$root)))
    } else {
      capture.output(invisible(try(uci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - alpha / 2,
                                                  interval = c(llim, ulim))$root)))
    }
    if(is.null(uci)){
      uci <- ulim
    }
  }

  if(is.nan(ltempPval)) {
    lci <- llim
  } else {
    lci <- NULL
    if(boot) {
      capture.output(invisible(try(lci <- uniroot(f = function(x) tibBootQuantile(x, theta, bootsamps, lower, upper, c, delta = delta) - (1 - alpha / 2),
                                                  interval = c(llim, ulim))$root)))
    } else {
      capture.output(invisible(try(lci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - (1 - alpha / 2),
                                                  interval = c(llim, ulim))$root)))
    }
    if(is.null(lci)) {
      lci <- llim
    }
  }

  ci <- c(lci, uci)
  return(ci)
}

ptruncNorm <- function(mu, x, sd, l, u) {
  ptruncnorm(x, l, u, mu, sd)
}

tibBootQuantile <- function(param, theta, samp, lower, upper, c, delta) {
  samp <- c * samp + param
  bootcdf <- (mean(theta <= samp & samp <= upper) + delta) / (mean(lower <= samp & samp <= upper) + delta)
  return(bootcdf)
}





