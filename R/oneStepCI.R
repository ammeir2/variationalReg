#' @export
oneStepCI <- function(object) {
  tCoef <- object$sampCoef
  n <- object$sampleSize
  p <- length(object$suffStat)
  selected <- object$selected
  ysig <- object$ysig
  cilevel <- object$cilevel

  condSamp <- rbind(as.numeric(object$suffStat), object$suffSamp)
  suffsamp <- condSamp / n
  forQuantiles <- apply(suffsamp, 2, function(x) x - mean(x))
  variance <- var(sqrt(n) * forQuantiles)
  A <- diag(p)
  A[selected, ] <- variance[selected, ]
  Ainv <- solve(A)
  forQuantiles <- forQuantiles %*% Ainv * ysig^2
  obs <- forQuantiles[1, ]
  forQuantiles <- forQuantiles[-1, ]
  quantiles <- apply(forQuantiles[, selected, drop = FALSE], 2,
                     function(x) quantile(x, probs = c(1 - cilevel / 2, cilevel / 2)))
  oneStep <- tCoef + obs[selected]
  weirdCI <- matrix(nrow = sum(selected), ncol = 2)
  for(i in 1:nrow(weirdCI)) {
    weirdCI[i, ] <- oneStep[i] - quantiles[, i]
  }
  return(weirdCI)
}
