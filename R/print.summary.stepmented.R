`print.summary.stepmented` <-
function(x, short = x$short, var.diff = x$var.diff, 
    digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"),...){
    cat("\n\t***Regression Model with Step Relationship(s)***\n\n")
    cat( "Call: \n" )
    print( x$call )
    cat("\nEstimated Jump-Point(s):\n ")
    est.psi<-x$psi[,c("Est.","St.Err"),drop=FALSE]
    #browser()
    est.psi[,1]<-x$psi.rounded[1,]

    rownames(est.psi)<-rownames(x$psi)
    print(round(est.psi,3)) #era "signif(,4)"
    if(!is.null(attr(x$psi.rounded,"break.dates"))) {
      cat("\ncorrisponding dates: ", attr(x$psi.rounded,"break.dates"),"\n")
      #cat(" ", attr(o$psi.rounded,"break.dates"),"\n")
    }
    
    if(short){ 
    cat("\nDifference-in-levels parameter(s):\n")
    #print(x$Ttable[(nrow(x$Ttable)-nrow(x$psi)+1):nrow(x$Ttable),])}
    nome<-rownames(x$psi)
    #nome<-as.character(parse("",text=nome))
    #aa<-grep("U",rownames(x$Ttable))
    #bb<-unlist(sapply(nome,function(xx){grep(xx,rownames(x$Ttable))},simplify=FALSE,USE.NAMES=FALSE))
    #cc<-intersect(aa,bb) #indices of diff-slope parameters
    nomiU<-rownames(x$gap)
    #idU<-match(nomiU,rownames(x$Ttable))
    print(x$Ttable[nomiU,])
      } else {cat("\nCoefficients of the linear terms:\n")
        if(is.null(dim(x$Ttable))){
        print(x$Ttable)
        #printCoefmat(matrix(x$Ttable,nrow=1,ncol=4,dimnames=list(" ",names(x$Ttable))),has.Pvalue=FALSE)
        } else {
        printCoefmat(x$Ttable, digits = digits, signif.stars = signif.stars,na.print = "NA", ...)
        }
        
        }
if("summary.lm"%in%class(x)){ #for lm
    if(var.diff){
    for(i in 1:length(x$sigma.new)){
    cat("\nResidual standard error ",i,":", format(signif(x$sigma.new[i], 
        digits)), "on", x$df.new[i], "degrees of freedom")}
    cat("\n")    
    } else {
    cat("\nResidual standard error:", format(signif(x$sigma, 
        digits)), "on", x$df[2], "degrees of freedom\n")}
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
        cat(",  Adjusted R-squared:", formatC(x$adj.r.squared, 
            digits = digits), "\n")}
        }
if("summary.glm"%in%class(x)){ #for glm
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n", apply(cbind(paste(format.default(c("Null", 
            "Residual"), width = 8, flag = ""), "deviance:"), 
            format(unlist(x[c("null.deviance", "deviance")]), 
                digits = max(5, digits + 1)), " on", format(unlist(x[c("df.null", 
                "df.residual")])), " degrees of freedom\n"), 
            1, paste, collapse = " "), "AIC: ", format(x$aic, 
            digits = max(4, digits + 1)), "\n", sep = "")
        }
if("summary.Arima"%in%class(x)){#for Arima 
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS") 
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), 
            ",  log likelihood = ", format(round(x$loglik, 2)), 
            ",  aic = ", format(round(x$aic, 2)), "\n", sep = "")
    else cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), 
        ",  part log likelihood = ", format(round(x$loglik, 2)), 
        "\n", sep = "")
}
    #browser()
invisible(x) 
#n.boot<-length(na.omit(..$psi.history$all.ss))
if(x$n.boot>0){
  cat("\nBoot restarting based on", x$n.boot, "samples.") #if(x$conv.warn) "*not*" else NULL , "attained in",x$it,"iter. (rel. change",paste(signif(x$epsilon,5),")\n",sep=""))
  cat("\nNumber of iterations in the last fit:",x$it, "(rel. change", paste(signif(x$epsilon,5),")\n",sep=""))
} else {
  cat("\nConvergence",if(x$conv.warn) "*not*" else NULL , "attained in",x$it,"iterations (rel. change",paste(signif(x$epsilon,5),")\n",sep=""))   
  }
}

