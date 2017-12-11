polyhedralMS <- function(X, y, ysig, selected, Eta = NULL, level = 0.95, computeCI = TRUE) {
  p <- ncol(X)
  suffStat <- t(X) %*% y
  sigma <- t(X) %*% X * ysig^2

  if(is.null(Eta)) {
    Eta <- matrix(0, ncol = sum(selected), nrow = ncol(X))
    Eta[selected, ] <- solve(t(X[, selected]) %*% X[, selected])
  }

  # Computing constraints -----------------------
  nconstraints <- 2 * sum(selected) * (p - sum(selected))
  alphaMat <- matrix(nrow = nconstraints, ncol = ncol(Eta))
  Ay <- numeric(nconstraints)
  etavars <- apply(Eta, 2, function(x) t(x) %*% sigma %*% x)
  whichSelected <- which(selected)
  notSelected <- which(!selected)

  constraintInd <- 0
  Arow <- rep(0, ncol(X))
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
    if(computeCI) polyCI[i, ] <- findPolyCIlimits(theta, etaSigma, Vminus, Vplus, 1 - level)
    pval[i] <- ptruncnorm(theta, Vminus, Vplus, 0, sqrt(etaSigma))
    pval[i] <- 2 * min(1 - pval[i], pval[i])
  }

  return(list(pval = pval, ci = polyCI))
}


findPolyCIlimits <- function(theta, etaSigma, lower, upper, alpha) {
  llim <- theta
  ltempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
  while(ltempPval < 1 - alpha / 2) {
    llim <- llim - sqrt(etaSigma) * 0.1
    ltempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
    if(is.nan(ltempPval)) {
      llim <- llim + sqrt(etaSigma) * 0.1
      break
    }
  }

  ulim <- theta
  utempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
  while(utempPval > alpha / 2) {
    ulim <- ulim + sqrt(etaSigma) * 0.1
    utempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
    if(is.nan(utempPval)) {
      ulim <- ulim - sqrt(etaSigma) * 0.1
      break
    }
  }

  if(is.nan(utempPval)) {
    uci <- ulim
  } else {
    uci <- NULL
    capture.output(invisible(try(uci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - alpha / 2,
                                                interval = c(llim, ulim))$root)))
    if(is.null(uci)){
      uci <- ulim
    }
  }

  if(is.nan(ltempPval)) {
    lci <- llim
  } else {
    lci <- NULL
    capture.output(invisible(try(lci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - (1 - alpha / 2),
                                                interval = c(llim, ulim))$root)))
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


