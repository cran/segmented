summary.segmented.lme<-function(object, .vcov=NULL, digits = max(3, getOption("digits") - 3), ...){
  #quale misura di varianza residua ? pi? approriata e deve restituire? da lme.fit o lme.fit.noG
  #stesso discorso vale per residui, logLik (and friends)
  #browser()
  x<- object$lme.fit
  dd <- x$dims
  LL<- as.numeric(logLik(object))
  
  #Quali residui?
  #resd <- resid(object$lme.fit, type = "pearson")
  resd <- resid(object$lme.fit.noG, type = "pearson")
  if (length(resd) > 5) {
    resd <- quantile(resd, na.rm = TRUE)
    names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
  }
  
  isLM<- inherits(x, "lm")
  
  fixeff <- if(isLM) coef else fixed.effects
  
  if(isLM){
    cat("Segmented model with varying breakpoint and slope difference\n")
  } else {
    cat("Segmented mixed-effects model fit by ")
    cat(if (x$method == "REML") "REML\n" else "maximum likelihood\n")
  }
  
  print(data.frame(AIC = AIC(object), BIC = BIC(object), logLik = LL, row.names = " "), ...)
  
  #      if(!is.null(object$history.boot.restart)) {
  #         n.sol<-length(unique(object$history.boot.restart[,"psi"]))
  #         cat("\nFit based on", nrow(object$history.boot.restart), "boot samples,", n.sol, "different solutions found\n")
  #         }
  if(!is.null(object$history.boot.restart)) {
    n.sol<-length(unique(object$history.boot.restart[,"psi"]))
    cat(" Bootstrap restarting on", nrow(object$history.boot.restart), "samples; ", n.sol, "different solution(s)\n")
  }
  
  cat("\n")
  
  if(!isLM){
    print(summary(x$modelStruct), sigma = x$sigma, reEstimates = x$coef$random,...) #  verbose = verbose, ...)
    cat("Fixed effects:\n")
  } else {
    cat("Estimates:\n")
  }
  
  
  #    fixF <- x$call$fixed
  #    if (inherits(fixF, "formula") || is.call(fixF)) {
  #        cat(deparse(x$call$fixed), "\n")
  #    }   else {
  #        cat(deparse(lapply(fixF, function(el) as.name(deparse(el)))), "\n")
  #    }
  
  #browser()
  
  nomeZ <- object$namesGZ$nameZ
  nomiZ <- c(paste(nomeZ,":",sep=""), paste(":",nomeZ,sep=""))
  nomiZ <- unlist(sapply(nomiZ, function(.x) grep(.x, names(fixeff(object$lme.fit)), value=TRUE)))
  nomiZ<-c(nomeZ, nomiZ)
  nomiDiffSlope <- c("U", object$namesGZ$nomiUx)
  nomiPsi <-c("G0", object$namesGZ$nomiG)
  nomiTutti <- names(fixeff(object$lme.fit))
  nomiLin <- setdiff(nomiTutti, c(nomiZ, nomiDiffSlope, nomiPsi))
  nominoG<-setdiff(nomiTutti, nomiPsi)
  if(isLM){
    a<-summary(x)$coefficients
    a[,"t value"]<- round(a[,"Estimate"]/a[,"Std. Error"], digits=3)
    a[,"Pr(>|t|)"]<-round(2*pt(-abs(a[,"t value"]), df=x$df.residual),digits=5)

  } else {
    a<-summary(x)$tTable #
    a[nominoG,"Value"] <- fixeff(object$lme.fit.noG)[nominoG]
    a[,"t-value"]<- round(a[,"Value"]/a[,"Std.Error"], digits=3)
    a[,"p-value"]<-round(2*pt(-abs(a[,"t-value"]), df=a[,"DF"]),digits=5)
  }
  
  
  id.leftS <- nomiZ %in% rownames(a)
  
  nomiOrd<- if(any(id.leftS)) c(nomiLin, nomiZ, nomiDiffSlope, nomiPsi ) else c(nomiLin, nomiDiffSlope, nomiPsi )
  a<-  a[nomiOrd,] 
  
  
  a["U", ncol(a)]<- NA #pValue per U
  
  a["G0",(ncol(a)-1):ncol(a)]<-NA
  if(!is.null(.vcov)){
    a[, "Std.Error"]<- sqrt(diag(.vcov))[nomiOrd]
    a[,"t-value"]<-round(a[,"Value"]/a[,"Std.Error"],digits=3)
    a[,"p-value"]<-round(2*pt(-abs(a[,"t-value"]), df=a[,"DF"]),4)
  }
  
  #a1<-rbind(a[nomiLin,],0,a[nomiZ,],0, a[nomiDiffSlope,],0, a[nomiPsi,])
  a1<-rbind(a[nomiLin,], NA, if(any(id.leftS)) a[nomiZ,] else NA, NA, a[nomiDiffSlope,], NA, a[nomiPsi,])
  #rownames(a1)<-c(nomiLin, "----", nomiZ, "----", nomiDiffSlope, "----", nomiPsi )
  rownames(a1)<-c(nomiLin, "-- leftS:", nomiZ, "-- diffS:", nomiDiffSlope, "-- break:", nomiPsi )
  if(!any(id.leftS)) a1[nomiZ, 1]<-0
  #a1[,ncol(a1)]<- round(a1[,ncol(a1)],4)
  #print(a1, na.print="", digits=4)
  #a1[,5]<- format.pval(a1[ ,5], eps = .001, digits = 2)
  print(a1, na.print="", digits=4, zero.print = "0")
  cat(" psi.link =", object$call$psi.link, "\n")
  #print(a1, zero.print="")
  
  cat("\nStandardized Within-Group Residuals:\n")
  print(resd, ...)
  
  #browser()
  
  
  cat("\nNumber of Observations:", x$dims[["N"]])
  cat("\nNumber of Groups: ")
  Ngrps <- dd$ngrps[1:dd$Q]
  if((lNgrps <- length(Ngrps)) == 1) {
    cat(Ngrps, "\n")
  } else {
    sNgrps <- 1:lNgrps
    aux <- rep(names(Ngrps), sNgrps)
    aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps, 
                                                   lNgrps))[!lower.tri(diag(lNgrps))])
    names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
    cat("\n")
    print(rev(Ngrps), ...)
  }
  
  #      if(!is.null(object$history.boot.restart)) {
  #         n.sol<-length(unique(object$history.boot.restart[,"psi"]))
  #         cat("\nFit based on", nrow(object$history.boot.restart), "boot samples,", n.sol, "different solutions found\n")
  #         }
  
  invisible(object)
  
}
