`print.segmented` <-
function(x,digits = max(3, getOption("digits") - 3),...){
#revisione 15/05/03; 24/02/04
if(is.null(x$psi)) x<-x[[length(x)]]
if(!"segmented"%in%class(x)) stop("a `segmented' object is requested")
cat( "Call: " )
print( x$call )
cat("\nMeaningful coefficients of the linear terms:\n")
#print(x$coef[(1:(length(x$coef)-length(x$psi[,2])))])
iV<- -match(x$nameUV[[2]],names(coef(x)))
#iV<- -grep("psi.",names(coef(x)))#indices all but V
print(x$coef[iV])
cat("\n")
cat("Estimated Break-Point(s)",dimnames(x$psi)[[1]],":",
    format(signif(x$psi[,2],digits)),"\n")
if("glm"%in%class(x)){    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
        x$df.residual, "Residual\n")
    cat("Null Deviance:    ", format(signif(x$null.deviance, 
        digits)), "\nResidual Deviance:", format(signif(x$deviance, 
        digits)), "     AIC:", format(signif(x$aic, digits)), "\n")
    }
if("Arima"%in%class(x)){
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS") 
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), 
            ":  log likelihood = ", format(round(x$loglik, 2)), 
            ",  aic = ", format(round(x$aic, 2)), "\n", sep = "")
    else cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), 
        ":  part log likelihood = ", format(round(x$loglik, 2)), 
        "\n", sep = "")
    }
    invisible(x)
}

