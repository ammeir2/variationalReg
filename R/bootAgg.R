waldSelection <- function(X, y, precision = NULL, level = 0.05) {
  if(is.null(precision)) {
    precision <- solve(sigma)
  }
  testStat <- t(y) %*% precision %*% y
  pval <- pchisq(testStat, df = length(y), lower.tail = FALSE)
  return(pval < level)
}

residualConditionalBootstrap <- function(X, y, ysig,
                                         cilevel = 0.05,
                                         thresholdLevel = cilevel^2 / length(y),
                                         selectionF = waldSelection, ...) {

}
