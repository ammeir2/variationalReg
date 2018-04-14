#' @export
univNormMLE <- function(y, sd, threshold) {
  if(length(threshold) == 1) {
    threshold <- abs(threshold)
    threshold <- c(-threshold, threshold)
  }

  fit <- nlm(univNormTruncDens, p = y, y = y, sd = sd, threshold = threshold)
  return(fit$estimate)
}

univNormTruncDens <- function(mu, y, sd, threshold) {
  dens <- dnorm(y, mu, sd, log = TRUE)
  prob <- pnorm(threshold[1], mu, sd, lower.tail = TRUE)
  prob <- prob + pnorm(threshold[2], mu, sd, lower.tail = FALSE)
  return(-dens + log(prob))
}
