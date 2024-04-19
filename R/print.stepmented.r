`print.stepmented` <- function(x, digits = max(3, getOption("digits") - 3),...){
#17/11/21!!! :))
        if(is.null(x$psi)) x<-x[[length(x)]]
        if(is.null(x$psi)) stop("the object does not include the 'psi' component")
        #if(!"segmented"%in%class(x)) stop("a `segmented' object is requested")
        cat( "Call: " )
        print( x$call )
        cat("\nCoefficients of the linear terms:\n")
        
        #browser()
        
        #iV<- -match(x$nameUV[[2]],names(coef(x)))#iV<- -grep("psi.",names(coef(x)))#indices all but V
        nomiPsi <- gsub("V", "psi", x$nameUV$V)
        coeff <- coef(x)
        iV<- -match(nomiPsi,names(coeff))#iV<- -grep("psi.",names(coef(x)))#indices all but V
        coeff <- if(any(is.na(coeff[iV]))) coeff else coeff[iV]  
        print.default(format(coeff, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
        cat("Estimated Jump-Point(s):\n")
        #a<-as.vector(x$psi[,"Est."]) 
        a<- as.vector(x$psi.rounded["inf [",])
        names(a)<-colnames(x$psi.rounded)
        print.default(a, digits+2, print.gap=2)
        if(!is.null(attr(x$psi.rounded,"break.dates"))) {
                cat("\ncorrisponding dates: ", attr(x$psi.rounded,"break.dates"),"\n")
                #cat(" ", attr(o$psi.rounded,"break.dates"),"\n")
        }
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

