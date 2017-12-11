summary.varMS <- function(object, ci_level = NULL, ci_types = NULL, verbose = TRUE, ...) {
  if(is.null(cilevel)) {
    cilevel <- object$cilevel
  }

  if(is.null(ci_types)) {
    ci_types <- c("naive", "naiveBoot")
    if(!is.null(object$varBootCI)) {
      ci_types <- c(ci_types, "varBoot")
    }
  }

  if("hybrid" %in% ci_types) {
    ci_types <- union(ci_types, "naive", "varBoot")
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

  if("naive" %in% ci_types) {
    naiveCI <- computeNaiveCI(object$suffStat, object$suffCov, object$selected, object$ysig, ci_level)
    cis[[slot]] <- naiveCI
    slot <- slot + 1
  }
}
