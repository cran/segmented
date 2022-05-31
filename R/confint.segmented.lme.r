confint.segmented.lme <- function(object, parm, level = 0.95, obj.boot, ...) {

  ci.boot <- function(m, conf.level = 0.95) {
    # computes three boot CI m: object returned by bootNP()
    est.orig <- m$coef[1, ]
    se.orig <- m$se[1, ]
    zalpha <- -qnorm((1 - conf.level)/2)
    EST <- m$coef[-1, ]
    # percentile
    CIt <- CIN <- CIperc <- apply(EST, 2, quantile, prob = c((1 - conf.level)/2,
                                                             (conf.level + (1 - conf.level)/2)), na.rm = TRUE)
    # Normal-based
    SE <- apply(EST, 2, sd, na.rm = TRUE)
    CIN[1, ] <- est.orig - zalpha * SE
    CIN[2, ] <- est.orig + zalpha * SE
    # t-boot
    Tdistr <- (EST - matrix(m$coef[1, ], ncol = length(est.orig), nrow = nrow(m$coef) -
                              1, byrow = TRUE))/m$se[-1, ]
    quantT <- apply(Tdistr, 2, quantile, prob = c((1 - conf.level)/2, (conf.level +
                                                                         (1 - conf.level)/2)), na.rm = TRUE)
    CIt[1, ] <- est.orig - quantT[2, ] * se.orig
    CIt[2, ] <- est.orig - quantT[1, ] * se.orig
    
    ris <- list(norm = CIN, perc = CIperc, t = CIt)
    ris
  }
  
  opz <- list(...)
  if (missing(obj.boot)) {
    if (is.null(opz$B)) {
      r <- object$lme.fit$varFix
      SE <- sqrt(diag(r))
      est <- object$lme.fit$coef$fixed
      zalpha <- -qnorm((1 - level)/2)
      CIN <- rbind(est - zalpha * SE, est + zalpha * SE)
      rownames(CIN) <- paste(100 * c((1 - level)/2, (level + (1 - level)/2)), "%", sep = "")
    } else {
      obj.boot <- vcov.segmented.lme(object, B = opz$B, ret.b=TRUE, seed = opz$seed, it.max.b = opz$it.max.b)
      CIN <- ci.boot(obj.boot, level)
    }
  } else {
    CIN <- ci.boot(obj.boot, level)
  }
  if(!missing(parm)){
    if(is.character(parm)) {
      parm <- match(parm, names(fixef(object)))
      if(any(is.na(parm))) stop("invalid names in parm")
    }
    if(is.list(CIN)) CIN <- lapply(CIN, function(x)x[,parm, drop=FALSE]) else CIN<-CIN[,parm,drop=FALSE]
  }
  return(CIN)
}


