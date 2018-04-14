plotEstimates <- function(estimates, CIs, offset = 0.12) {
  p <- length(estimates[[1]])
  nests <- length(estimates)
  ncis <- length(CIs)
  plotList <- list()
  for(i in 1:max(nests, ncis)) {
    if(nests >= i) {
      est <- estimates[[i]]
      estname <- names(estimates)[i]
    } else {
      est <- rep(NA, p)
      estname <- NA
    }

    if(ncis >= i) {
      ci <- CIs[[i]]
      ciname <- names(CIs)[i]
    } else {
      ci <- matrix(NA, nrow = p, ncol = 2)
      ciname <- NA
    }

    plotList[[i]] <- data.frame(variable = 1:p, estimate = est,
                                lci = ci[, 1], uci = ci[, 2],
                                estimate_type = estname, ci_type = ciname,
                                offset = i)
  }

  offsetMulti <- offset
  plotdat <- do.call("rbind", plotList)

  # plotting
  plotdat$variable <- plotdat$variable + plotdat$offset * offsetMulti
  ggplot(subset(plotdat, !is.na(ci_type))) +
    geom_segment(aes(x = variable, xend = variable,
                     y = lci, yend = uci, col = ci_type, linetype = ci_type)) +
    theme_bw() + geom_hline(yintercept = 0) +
    geom_point(data = subset(plotdat, !is.na(estimate_type)),
               aes(x = variable, y = estimate, shape = estimate_type)) +
    scale_color_discrete(name = "CI Type") +
    # scale_color_discrete(name = "Estimate Type") +
    scale_shape_discrete(name = "Estimate Type") +
    scale_linetype_discrete(name = "CI Type") +
    ylab("Estimates / CIs") + xlab("Variable")
}
