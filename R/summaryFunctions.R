summary.varMS <- function(object, ci_level = NULL, ci_types = NULL, verbose = TRUE,
                          naiveBound = TRUE, ...) {
  if(is.null(ci_level)) {
    cilevel <- object$cilevel
  } else {
    cilevel <- ci_level
  }

  if(is.null(ci_types)) {
    ci_types <- c("naive", "naiveBoot")
    if(!is.null(object$varBootCI)) {
      ci_types <- c(ci_types, "varBoot")
    }
  }

  cis <- vector('list', length(ci_types))
  slot <- 1

  if("varBoot" %in% ci_types & is.null(object$varBoot)) {
    XtX <- object$suffCov / object$ysig^2
    XtXinv <- solve(XtX)
    suffPrecision <- XtXinv / object$ysig^2
    varBoot <- computeVarBootSample(object$suffSamp, suffPrecision, object$suffCov, XtXinv, object$threshold,
                                    verbose = verbose)
    object$varBoot <- varBoot
  }

  naiveCI <- computeNaiveCI(as.numeric(object$suffStat), object$suffCov, as.numeric(object$selected), object$ysig, cilevel)
  if("naive" %in% ci_types) {
    cis[[slot]] <- naiveCI
    slot <- slot + 1
  }

  if("naiveBoot" %in% ci_types) {
    thresholdCoef <- object$naive
    thresholdCoef[object$sampCoef == 0] <- 0
    naiveBootCI <- computeBootCI(object$naiveBoot, thresholdCoef, object$naive, cilevel)
    if(naiveBound) {
      naiveBootCI[, 1] <- pmin(naiveBootCI[, 1], naiveCI[, 1])
      naiveBootCI[, 2] <- pmax(naiveBootCI[, 2], naiveCI[, 2])
    }
    cis[[slot]] <- naiveBootCI
    slot <- slot + 1
  }

  if("varBoot" %in% ci_types) {
    varBootCI <- computeBootCI(object$varBoot, object$sampCoef, object$mEst, cilevel)
    if(naiveBound) {
      varBootCI[, 1] <- pmin(varBootCI[, 1], naiveCI[, 1])
      varBootCI[, 2] <- pmax(varBootCI[, 2], naiveCI[, 2])
    }
    cis[[slot]] <- varBootCI
    slot <- slot + 1
  }

  if("poly" %in% ci_types) {
    polyfit <- polyhedralMS(object$suffStat, object$suffCov, object$ysig,
                            object$selected, Eta = NULL,
                            level = 1 - cilevel)
    cis[[slot]] <- polyfit$ci
  }

  return(cis)
}
