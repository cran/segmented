
#--generic function
segmented<-function(obj, Z, psi, W, it.max=20, toll=0.0001, visual=FALSE, last=TRUE, ...){
            UseMethod("segmented")
            }


#--default method
segmented.default<-function(obj, Z, psi, W, it.max=20, toll=0.0001, visual=FALSE, last=TRUE, ...){
            stop("No default method for segmented")
            }

#--method for lm objects
segmented.lm <-
#revisione 13/05/03; 22/09/03; 02/10/03; 9/10/03; 3/11/03; 27/11/03; 20/01/04
function(obj, Z, psi, W, it.max=20, toll=0.0001, visual=FALSE, last=TRUE,...){
if(is.data.frame(eval(obj$call$data))){
    attach(eval(obj$call$data))
    on.exit(detach(eval(obj$call$data)))}
    if(is.null(obj$y) || is.null(dim(obj$x))) obj<-update(obj,x=TRUE,y=TRUE) 
    y<-obj$y
    if(is.matrix(Z)) name.Z<-colnames(Z)
        else name.Z<-deparse(substitute(Z))
#Da aggiungere da qua....
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
        obj<-update(obj,.~.+U)
        obj$psi<-psi
        return(obj)}
    initial<-psi
    it<-1
    epsilon<-10
    XREG<-obj$x    
    o<-if(!is.null(obj$offset)) obj$offset else NULL
    w<-if(is.null(obj$weights)) rep(1,length(y)) else obj$weights 
    obj0<-obj
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
                    if(it==1) cat(0," ",formatC(dev.new,3,format="f"),"","(No breakpoint(s))","\n")
                              spp<-if(it<10) "" else NULL
                              cat(it,spp,"",formatC(dev.new,3,format="f"),"\n")}
            epsilon<- (dev.new-dev.old)/dev.old
            obj$epsilon<-epsilon
            it<- it+1
            obj$it<-it
            class(obj)<-c("segmented",class(obj))
            list.obj[[length(list.obj)+ifelse(last==TRUE,0,1)]]<-obj
            if(it>it.max) break
            }
if(it>it.max) warning("max number of iterations attained",call.=FALSE)
Vxb <-t(t(V)*beta.c)
X[,(dim(X)[2]+1-dim(Vxb)[2]):dim(X)[2]]<-Vxb
rownames(X)<-NULL
#
nameVxb<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="")
colnames(X)[(ncol(X)+1-ncol(Vxb)):ncol(X)]<-nameVxb
nameU<- if(k==1) paste("U",".",name.Z,sep="") else paste("U", 1:k, ".", name.Z, sep="")
colnames(X)[(ncol(X)-2*ncol(Vxb)+1):(ncol(X)-ncol(Vxb))]<-nameU
obj0$model[,nameU]<-U
obj0$model[,nameVxb]<-Vxb
obj.final<- update(obj0,formula=
as.formula(
paste(deparse(formula(obj0)),paste("`",nameU,"`",collapse="+",sep=""),paste("`",nameVxb,"`",collapse="+",sep=""),sep="+"))
,data=obj0$model,...)
Cov<-summary(obj.final)$cov.unscaled*summary(obj.final)$sigma^2
testo<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="")
if(length(grep(":",name.Z))>0) testo<-paste("`",testo,"`",sep="")
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
#revisione 22/09/03; 02/10/03; 03/11/03; 27/11/03; 20/01/04
function(obj, Z, psi, W, it.max=20, toll=0.0001, visual=FALSE, last=TRUE, ...){
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
    obj0<-obj
    list.obj<-list(obj)
        while(abs(epsilon)>toll){
            U<-pmax((Z -PSI), 0) 
            V<-ifelse((Z >PSI), -1, 0)
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
                    if(it==1) cat(0," ",formatC(dev.new,3,format="f"),"","(No breakpoint(s))","\n")
                              spp<-if(it<10) "" else NULL 
                              cat(it,spp,"",formatC(dev.new,3,format="f"),"\n")}
            epsilon<- (dev.new-dev.old)/dev.old
            obj$epsilon<-epsilon
            it<- it+1
            obj$it<-it
            class(obj)<-c("segmented",class(obj))
            list.obj[[length(list.obj)+ifelse(last==TRUE,0,1)]]<-obj
            if(it>it.max) break
            }
if(it>it.max) warning("max number of iterations attained",call.=FALSE)
Vxb <-t(t(V)*beta.c)
X[,(dim(X)[2]+1-dim(Vxb)[2]):dim(X)[2]]<-Vxb
rownames(X)<-NULL
nameVxb<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="")#nuovo
colnames(X)[(ncol(X)+1-ncol(Vxb)):ncol(X)]<-nameVxb
nameU<- if(k==1) paste("U",".",name.Z,sep="") else paste("U", 1:k, ".", name.Z, sep="") #aggiunto
colnames(X)[(ncol(X)-2*ncol(Vxb)+1):(ncol(X)-ncol(Vxb))]<-nameU #aggiunto
obj0$model[,nameU]<-U
obj0$model[,nameVxb]<-Vxb
obj.final<- update(obj0,formula=
as.formula(
paste(deparse(formula(obj0)),paste("`",nameU,"`",collapse="+",sep=""),paste("`",nameVxb,"`",collapse="+",sep=""),sep="+"))
,data=obj0$model,...)
Cov<-summary(obj.final)$cov.scaled
testo<- if(k==1) paste("psi",".",name.Z,sep="") else paste("psi", 1:k, ".", name.Z, sep="") #nuovo
if(length(grep(":",name.Z))>0) testo<-paste("`",testo,"`",sep="")
vv<- if(k==1) Cov[testo,testo] else diag(Cov[testo,testo])#nuovo
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

summary.segmented<-function(object, short=FALSE, ...){
#revisione 13/05/03;7/10/03;28/11/03
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
        summ$short<-short #AGG
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
        summ$short<-short #AGG
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
        object$short<-short #AGG
        summ<-object 
        class(summ) <- "summary.segmented"
        return(summ)}
}

#--------------

print.summary.segmented<-function(x,digits = max(3, getOption("digits") - 3),short=x$short,...){
#revisione 15/05/03;7/10/03;28/11/03
# con segmented.lm()
    cat("\n\t***Regression Model with Segmented Relationship(s)***\n\n")
    cat( "Call: \n" )
    print( x$call )
    cat("\nEstimated Break-Point(s):\n ")
    print(signif(x$psi[,-1],4))
    cat("\nt value for the gap-variable(s) V: ",x$gap[,3],"\n")
if(any(abs(x$gap[,3])>1.96)) cat("    Warning: some coefficient of the gap-variable is significant at 0.05 level\n")
    if(short){ 
    cat("\nDifference-in-slopes parameter(s):\n")
    print(x$Ttable[(nrow(x$Ttable)-nrow(x$psi)+1):nrow(x$Ttable),])}
    else {cat("\nMeaningful coefficients of the linear terms:\n")
        print(x$Ttable)}
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

#-------------------------

slope.segmented<-function(ogg, level=0.95){
#2/10/03; 21/10/03; 4/11/03
#Ad ogni chiamata di grep() è stato aggiunto l'argomento extended=F per consentire Z=log(x)
#returns Est., St.Err., t value and Conf Interv (1-level) for the slopes
#      of the variables having a segmented relationship in the model.
#ogg: a segmented object returned by segmented().
        if(!"segmented"%in%class(ogg)) stop("A segmented model is needed")
        nome<-rownames(ogg$psi)
        NOME<-names(coef(eval(ogg$call$obj)))
        index<-NULL
        nome1<-as.character(parse("",text=nome))
        cc<-grep(nome1, NOME, extended=FALSE) 
        if(length(nome1)==1) { 
                index<-grep(nome1, names(coef(ogg)), extended=FALSE)
                index<-list(index[1:(length(index)-nrow(ogg$psi))])
                    
                    if(length(cc)>1) {
                        index<-unlist(index)
                        index<-cbind(index[-length(index)],index[length(index)])
                        index<-apply(index,1,list)}
                            if(length(cc)==1 && !is.null(ogg$call$W)) { #quando c'è W e cc=1
                                index<-unlist(index)
                                index<-cbind(index[1],index[-1])
                                index<-apply(index,1,list)
                                        }
                                }
            else {
                for(i in 1:length(nome1)){
                        aa<-grep(nome1[i], names(coef(ogg)), extended=FALSE)
                        index[[i]]<-list(aa[1:(length(aa)-1)])
                        }
                }
        Ris<-list()   
        for(i in 1:length(index)){
            ind<-unlist(index[[i]])
            M<-matrix(1,length(ind),length(ind))
            M[row(M)<col(M)]<-0
            cof<-ogg$coef[ind]
            covv<- if("Arima"%in%class(ogg)){
                    ogg$var.coef[ind,ind] } else 
                    {if("glm"%in%class(ogg)) summary.glm(ogg)$cov.scaled[ind,ind]
                        else summary.glm(ogg)$cov.unscaled[ind,ind]*summary(ogg)$sigma^2
                        }
            cof.out<-M%*%cof 
            cov.out<-M%*%covv%*%t(M)
            se.out<-sqrt(diag(cov.out))
            k<-abs(qnorm((1-level)/2))*se.out
            ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
            cin<-paste("CI","(",level*100,"%",")",c(".l",".u"),sep="")
            dimnames(ris)<-list(paste("slope", 1:length(ind), sep=""),c("Est.","St.Err.","t value",cin[1],cin[2]))
            Ris[[nome1[i]]]<-ris}
            if(any(is.na(names(Ris)))){
                nn1<-paste(":",deparse(ogg$call$Z),sep="")
                nn2<-paste(deparse(ogg$call$Z),":",sep="")
                l<-list(grep(nn1, NOME, extended=FALSE), grep(nn2, NOME, extended=FALSE))
                f<-sapply(l,function(x){length(x)!=0})
#                if(!is.null(ogg$call$W) || length(l[f])>0) {
                nomeW<-if(!is.null(ogg$call$W)) 
                        paste(deparse(ogg$call$W),1:length(index),":",deparse(ogg$call$Z),sep="")
                            else names(ogg$coef)[unlist(l[f])]
            names(Ris)<-nomeW #}
                    }
            Ris
            }
