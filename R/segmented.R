
#--generic function
segmented<-function(obj, Z, psi, W, it.max, toll, visual, last){
            UseMethod("segmented")
            }


#--default method
segmented.default<-function(obj, Z, psi, W, it.max, toll, visual, last){
            stop("No default method for segmented")
            }

#--method for lm objects
segmented.lm <-
#revisione 13/05/03; 22/09/03; 02/10/03
function(obj, Z, psi, W, it.max=20, toll=0.0001, visual=FALSE, last=TRUE){
if(is.data.frame(eval(obj$call$data))){
    attach(eval(obj$call$data))
    on.exit(detach(eval(obj$call$data)))}
    if(is.null(obj$y) || is.null(dim(obj$x))) obj<-update(obj,x=TRUE,y=TRUE) 
    y<-obj$y
    if(is.matrix(Z)) name.Z<-dimnames(Z)[[2]]
        else name.Z<-deparse(substitute(Z))
    if(!missing(W) && !is.matrix(Z)) {#L groups
                if(length(table(W))!=length(psi)) stop("number of levels of W and length(psi) do not match")
                        Z<-model.matrix(~ factor(W)-1)*Z
                        sel<-grep(name.Z,colnames(obj$x))
                        name.Z<-colnames(obj$x)[sel]}
    if(missing(W) && !is.matrix(Z)){#multi
                        Z<-matrix(rep(Z,length(psi)),ncol=length(psi))}
    k<-ncol(Z)
    dimnames(Z)<-list(NULL, rep("",k))
            if(ncol(Z)!=length(psi)) stop("number of Z variables and length(psi) don't match")
            if(nrow(Z)!=length(obj$y)) stop("length(Z) is different from length(obj$y)")
    PSI<- matrix(rep(psi, rep(nrow(Z),ncol(Z))), ncol=ncol(Z))
    if(it.max==0){
        U<-pmax((Z -PSI), 0)
        rownames(U)<-NULL
        obj<-update(obj,.~.+U)
        obj$psi<-psi
        return(obj)}
    initial<-psi
    it<-1
    epsilon<-10
    XREG<-obj$x    
    o<-if(!is.null(obj$offset)) obj$offset else NULL
    w<-if(is.null(obj$weights)) rep(1,length(y)) else obj$weights 
    list.obj<-list(obj)
        while(abs(epsilon)>toll){
            U<-pmax((Z -PSI), 0) 
            V<-ifelse((Z >PSI), -1, 0)
            dev.old<-sum(obj$residuals^2)
            X<-cbind(XREG,U,V)
            rownames(X)<-NULL
            if(ncol(V)==1) colnames(X)[(ncol(XREG)+1):ncol(X)]<-c("U","V") 
                else colnames(X)[(ncol(XREG)+1):ncol(X)]<-c(paste("U", 1:k, sep=""),paste("V", 1:k, sep=""))
            obj<-lm.wfit(x=X,y=y,w=w,offset=o)        
            if(k==1) beta.c<-coef(obj)["U"] else beta.c<-coef(obj)[paste("U", 1:k, sep="")]
            if(k==1) gamma.c<-coef(obj)["V"] else gamma.c<-coef(obj)[paste("V", 1:k, sep="")]
            psi.old<- psi
            psi<-psi.old+gamma.c/beta.c
            PSI<- matrix(rep(psi, rep(nrow(Z),ncol(Z))), ncol=ncol(Z))
            a<-apply((Z<PSI),2,all)
            b<-apply((Z>PSI),2,all)
            if(sum(a+b)!=0 || is.na(sum(a+b)))
                stop("(Some) estimated psi out of its range")
            obj$psi<-psi
            dev.new<-sum(obj$residuals^2)
                if(visual) {
                    if(it==1) cat(0,"",round(dev.old,3),"","(No breakpoint(s))","\n")
                              cat(it,"",round(dev.new,3),"\n")}
            epsilon<- (dev.new-dev.old)/dev.old
            obj$epsilon<-epsilon
            it<- it+1
            obj$it<-it
            class(obj)<-c("segmented",class(obj))
            list.obj[[length(list.obj)+ifelse(last==TRUE,0,1)]]<-obj
            if(it>it.max) break
            }
if(it>it.max) warning("Convergence attained with max iterations")
Vxb <-t(t(V)*beta.c)
X[,(dim(X)[2]+1-dim(Vxb)[2]):dim(X)[2]]<-Vxb
rownames(X)<-NULL
nameVxb<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="")
colnames(X)[(ncol(X)+1-ncol(Vxb)):ncol(X)]<-nameVxb
nameU<- if(k==1) paste("U",".",name.Z,sep="") else paste("U", 1:k, ".", name.Z, sep="")
colnames(X)[(ncol(X)-2*ncol(Vxb)+1):(ncol(X)-ncol(Vxb))]<-nameU
obj.final<- lm(y~X-1, offset=o,weights=w)
names(obj.final$coefficients)<- colnames(X)
Cov<-summary(obj.final)$cov.unscaled*summary(obj.final)$sigma^2
dimnames(Cov)<-list(colnames(X),colnames(X))
testo<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="")
vv<- if(k==1) Cov[testo,testo] else diag(Cov[testo,testo])
obj.final$psi<-cbind(initial, psi, sqrt(vv))
dimnames(obj.final$psi)[[1]] <-if(length(name.Z)==1) c(name.Z,rep("",k-1)) else name.Z
dimnames(obj.final$psi)[[2]] <-c("Initial","Est.","St.Err")
obj.final$it<-(it-1)
obj.final$epsilon<-epsilon
obj.final$call<-match.call()
class(obj.final)<- c("segmented","lm")
list.obj[[length(list.obj)+1]]<-obj.final
class(list.obj)<-"segmented"
if(last) list.obj<-list.obj[[length(list.obj)]]
return(list.obj)
}





#--method for glm objects
segmented.glm <-
#revisione 22/09/03; 02/10/03
function(obj, Z, psi, W, it.max=20, toll=0.0001, visual=FALSE, last=TRUE){
if(is.data.frame(obj$data)){
        attach(obj$data)
        on.exit(detach(obj$data))}
    if(is.null(obj$y) || is.null(dim(obj$x))) obj<-update(obj,x=TRUE,y=TRUE) 
    y<-obj$y
    if(is.matrix(Z)) name.Z<-dimnames(Z)[[2]]
        else name.Z<-deparse(substitute(Z))
    if(!missing(W) && !is.matrix(Z)) {#L groups
                if(length(table(W))!=length(psi)) stop("number of levels of W and length(psi) do not match")
                        Z<-model.matrix(~ factor(W)-1)*Z
                        sel<-grep(name.Z,colnames(obj$x))
                        name.Z<-colnames(obj$x)[sel]}
    if(missing(W) && !is.matrix(Z)){#multi
                        Z<-matrix(rep(Z,length(psi)),ncol=length(psi))}
    k<-ncol(Z)
    dimnames(Z)<-list(NULL, rep("",k))
            if(ncol(Z)!=length(psi)) stop("number of Z variables and length(psi) do not match")
            if(nrow(Z)!=length(obj$y)) stop("length(Z) is different from length(obj$y)")
    PSI<- matrix(rep(psi, rep(nrow(Z),ncol(Z))), ncol=ncol(Z))
    if(it.max==0){
        U<-pmax((Z -PSI), 0)
        rownames(U)<-NULL
        obj$model$U<-U
        obj<-update(obj,.~.+U,data=obj$model)
        obj$psi<-psi
        return(obj)}
    initial<-psi
    it<-1
    epsilon<-10
    XREG<-update(obj,x=TRUE)$x
    y<-obj$y
    fam<-family(obj)
    o<-obj$offset
    contr<-obj$control
    w<-obj$prior.weights
    list.obj<-list(obj)
        while(abs(epsilon)>toll){
            U<-pmax((Z -PSI), 0) 
            V<-ifelse((Z >PSI), -1, 0)
            #dimnames(V)<-list(NULL, 1:k) #forse superfluo
            #dimnames(U)<-list(NULL, 1:k) #forse superfluo
            dev.old<-obj$dev
            X<-cbind(XREG,U,V)
            rownames(X)<-NULL
            if(ncol(V)==1) colnames(X)[(ncol(XREG)+1):ncol(X)]<-c("U","V") 
                else colnames(X)[(ncol(XREG)+1):ncol(X)]<-c(paste("U", 1:k, sep=""),paste("V", 1:k, sep=""))
            obj<-glm.fit(x=X,y=y,offset=o,weights=w,family=fam,control=contr)
            if(k==1) beta.c<-coef(obj)["U"] else beta.c<-coef(obj)[paste("U", 1:k, sep="")]
            if(k==1) gamma.c<-coef(obj)["V"] else gamma.c<-coef(obj)[paste("V", 1:k, sep="")]
            psi.old<- psi
            psi<-psi.old+gamma.c/beta.c
            PSI<- matrix(rep(psi, rep(nrow(Z),ncol(Z))), ncol=ncol(Z))
            a<-apply((Z<PSI),2,all)
            b<-apply((Z>PSI),2,all)
            if(sum(a+b)!=0 || is.na(sum(a+b)))
                stop("(Some) estimated psi out of its range")
            obj$psi<-psi
            dev.new<-obj$dev 
                if(visual) {
                    if(it==1) cat(0,"",round(dev.old,3),"","(No breakpoint(s))","\n")
                              cat(it,"",round(dev.new,3),"\n")}
            epsilon<- (dev.new-dev.old)/dev.old
            obj$epsilon<-epsilon
            it<- it+1
            obj$it<-it
            class(obj)<-c("segmented",class(obj))
            list.obj[[length(list.obj)+ifelse(last==TRUE,0,1)]]<-obj
            if(it>it.max) break
            }
if(it>it.max) warning("Convergence attained with max iterations")
Vxb <-t(t(V)*beta.c)
X[,(dim(X)[2]+1-dim(Vxb)[2]):dim(X)[2]]<-Vxb
rownames(X)<-NULL
nameVxb<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="")
colnames(X)[(ncol(X)+1-ncol(Vxb)):ncol(X)]<-nameVxb
nameU<- if(k==1) paste("U",".",name.Z,sep="") else paste("U", 1:k, ".", name.Z, sep="") 
colnames(X)[(ncol(X)-2*ncol(Vxb)+1):(ncol(X)-ncol(Vxb))]<-nameU 
obj.final<- glm(y~X-1, offset=o,weights=w,family=fam,control=contr)
names(obj.final$coefficients)<- colnames(X)
Cov<-summary(obj.final)$cov.scaled
dimnames(Cov)<-list(colnames(X),colnames(X))
testo<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="") 
vv<- if(k==1) Cov[testo,testo] else diag(Cov[testo,testo])
obj.final$psi<-cbind(initial, psi, sqrt(vv))
dimnames(obj.final$psi)[[1]] <-if(length(name.Z)==1) c(name.Z,rep("",k-1)) else name.Z
dimnames(obj.final$psi)[[2]] <-c("Initial","Est","St.Err")
obj.final$it<-it-1
obj.final$epsilon<-epsilon
obj.final$call<-match.call()
class(obj.final)<- c("segmented","glm","lm")
list.obj[[length(list.obj)+1]]<-obj.final
class(list.obj)<-"segmented"
if(last) list.obj<-list.obj[[length(list.obj)]]
return(list.obj)
}


##Funzioni metodo

print.segmented<-function(x,digits = max(3, getOption("digits") - 3),...){
#revisione 15/05/03
#ok con segmented.lm()
if(is.null(x$psi)) x<-x[[length(x)]]
if(!"segmented"%in%class(x)) stop("a `segmented' object is requested")
cat( "Call: " )
print( x$call )
cat("\nMeaningful coefficients of the linear terms:\n")
print(x$coef[(1:(length(x$coef)-length(x$psi[,2])))])
cat("\n")
cat("Estimated Break-Point(s) for the variable(s)",dimnames(x$psi)[[1]],":",
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

summary.segmented<-function(object, ...){
#revisione 13/05/03;7/10/03
# con segmented.lm()?
if(is.null(object$psi)) object<-object[[length(object)]]
iV<-(1:(length(object$coef)-length(object$psi[,2])))#indices of all but the Vs
beta.c<- if(length(object$psi[,1])==1) object$coef["U"] else object$coef[paste("U", 1:length(object$psi[,1]), sep="")]
    if("lm"%in%class(object) && !"glm"%in%class(object)){
        summ <- c(summary.lm(object, ...), object["psi"])
        summ$Ttable<-summ$coefficients[iV,-4]
        coeff<-summ$coefficients[,1]
        v<-summ$coefficients[,2]
        summ$gap<-cbind(coeff[-iV]*beta.c,v[-iV]*beta.c,coeff[-iV]/v[-iV])
        dimnames(summ$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        summ[c("it","epsilon")]<-object[c("it","epsilon")]
        class(summ) <- c("summary.segmented", "summary.lm")
        return(summ)}
    if("glm"%in%class(object)){
        summ <- c(summary.glm(object, ...), object["psi"])
        summ$Ttable<-summ$coefficients[iV,-4]
        coeff<-summ$coefficients[,1]
        v<-summ$coefficients[,2]
        summ$gap<-cbind(coeff[-iV]*beta.c,v[-iV]*beta.c,coeff[-iV]/v[-iV])
        dimnames(summ$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        summ[c("it","epsilon")]<-object[c("it","epsilon")]
        class(summ) <- c("summary.segmented", "summary.glm")
        return(summ)}
    if("Arima"%in%class(object)){
        coeff<-object$coef
        v<-sqrt(diag(object$var.coef))
        Ttable<-cbind(coeff[iV],v[iV],coeff[iV]/v[iV])
        object$gap<-cbind(coeff[-iV]*beta.c,v[-iV]*beta.c,coeff[-iV]/v[-iV])
        dimnames(object$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        colnames(Ttable)<-c("Estimate","Std. Error","t value")
        object$Ttable<-Ttable
        summ<-object 
        class(summ) <- "summary.segmented"
        return(summ)}
}

#--------------

print.summary.segmented<-function(x,digits = max(3, getOption("digits") - 3),...){
#revisione 15/05/03;7/10/03
# con segmented.lm()
    cat("\n\t***Regression Model with Segmented Relationship(s)***\n\n")
    cat( "Call: \n" )
    print( x$call )
    cat("\nEstimated Break-Point(s):\n ")
    print(signif(x$psi[,-1],4))
    cat("\nt value for the gap-variable(s) V: ",x$gap[,3],"\n")
if(any(abs(x$gap[,3])>1.96)) cat("    Warning: some coefficient of the gap-variable is significant at 0.05 level\n")
    cat("\nMeaningful coefficients of the linear terms:\n")
        print(x$Ttable)
    cat("\n")
if("summary.lm"%in%class(x)){ #for lm
    cat("\nResidual standard error:", format(signif(x$sigma, 
        digits)), "on", x$df[2], "degrees of freedom\n")
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
        cat(",  Adjusted R-squared:", formatC(x$adj.r.squared, 
            digits = digits), "\n")}
        }

if("summary.glm"%in%class(x)){ #for glm
    cat("(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n", apply(cbind(paste(format.char(c("Null", 
            "Residual"), width = 8, flag = ""), "deviance:"), 
            format(unlist(x[c("null.deviance", "deviance")]), 
                digits = max(5, digits + 1)), " on", format(unlist(x[c("df.null", 
                "df.residual")])), " degrees of freedom\n"), 
            1, paste, collapse = " "), "AIC: ", format(x$aic, 
            digits = max(4, digits + 1)), "\n", sep = "")
        }

if(!"summary.lm"%in%class(x) && !"summary.glm"%in%class(x)){#for Arima 
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS") 
        cat("sigma^2 estimated as ", format(x$sigma2, digits = digits), 
            ",  log likelihood = ", format(round(x$loglik, 2)), 
            ",  aic = ", format(round(x$aic, 2)), "\n", sep = "")
    else cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), 
        ",  part log likelihood = ", format(round(x$loglik, 2)), 
        "\n", sep = "")
    }
invisible(x) 
cat("\nConvergence attained in",x$it,"iterations with relative change",x$epsilon,"\n")
}
