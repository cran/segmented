print.segmented.lme<- function (x, digits = max(3, getOption("digits") - 3), ...) {
  #datacall<- eval(x$call$obj)$call$data
  
  #datacall<- if(is.call(eval(x$call$obj))) eval(x$call$obj)$data else  eval(x$call$obj)$call$data
  datacall<- x$misc$datacall
  xx<-x
  LL<-x$lme.fit.noG$logLik
  x<-x$lme.fit 
  dd <- x$dims
  cat("Segmented linear mixed-effects model fit by ")
  cat(if (x$method == "REML") 
    "REML\n"
    else "maximum likelihood\n")
  cat("  Data:", datacall, "\n")
  if (!is.null(x$call$subset)) {
    cat("  Subset:", deparse(asOneSidedFormula(x$call$subset)[[2L]]), 
        "\n")
  }
  cat("  Log-", if (x$method == "REML") 
    "restricted-"
    else "", "likelihood (approx): ", format(LL), "\n", sep = "")
  if(!is.null(xx$history.boot.restart)) {
    n.sol<-length(unique(xx$history.boot.restart[,"psi"]))
    cat("  Bootstrap restarting on", nrow(xx$history.boot.restart), "samples;", n.sol, "different solution(s)\n")
  }
  #cat(" \n psi.link =", xx$call$psi.link, "\n")
  cat("\n")
  fixF <- x$call$fixed
  cat("Fixed:", deparse(if (inherits(fixF, "formula") || 
                            is.call(fixF) || is.name(fixF)) 
    x$call$fixed
    else lapply(fixF, function(el) as.name(deparse(el)))), "\n")
  print(fixef(xx), ...) #<-xx e' l'oggetto segmented.lme
  cat("  psi.link =", xx$call$psi.link)
  cat("\n\n")
  print(summary(x$modelStruct), sigma = x$sigma, ...)
  cat("Number of Observations:", dd[["N"]])
  cat("\nNumber of Groups: ")
  Ngrps <- dd$ngrps[1:dd$Q]
  if ((lNgrps <- length(Ngrps)) == 1) {
    cat(Ngrps, "\n")
  }
  else {
    sNgrps <- 1:lNgrps
    aux <- rep(names(Ngrps), sNgrps)
    aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps, 
                                                   lNgrps))[!lower.tri(diag(lNgrps))])
    names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
    cat("\n")
    print(rev(Ngrps), ...)
  }
  invisible(x)
}
