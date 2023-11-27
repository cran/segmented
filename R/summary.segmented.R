`summary.segmented` <-
function(object, short=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, ...){
    if(is.null(object$psi)) object<-object[[length(object)]]
#i seguenti per calcolare aa,bb,cc funzionano per lm e glm, da verificare con arima....
#    nome<-rownames(object$psi)
#    nome<-as.character(parse("",text=nome))
#    aa<-grep("U",names(coef(object)[!is.na(coef(object))]))
#    bb<-unlist(sapply(nome,function(x){grep(x,names(coef(object)[!is.na(coef(object))]))},simplify=FALSE,USE.NAMES=FALSE))
#    cc<-intersect(aa,bb) #indices of diff-slope parameters
#    iV<- -grep("psi.",names(coef(object)[!is.na(coef(object))]))#indices of all but the Vs
    if(!is.null(.vcov)) var.diff<-FALSE
    if(var.diff && length(object$nameUV$Z)>1) {
      var.diff<-FALSE
      warning(" 'var.diff' set to FALSE with multiple segmented variables", call.=FALSE)
    }
    nomiU <- object$nameUV$U
    nomiV <- object$nameUV$V
    idU<-match(nomiU,names(coef(object)[!is.na(coef(object))]))
    idV<-match(nomiV,names(coef(object)[!is.na(coef(object))]))
    beta.c<- coef(object)[nomiU]
    #per metodo default.. ma serve????
    #browser()
    if(all(is.na(object[["psi"]][,"St.Err"]))) {
      if(inherits(object, "lm")){
        R <- chol2inv(object$qr$qr)
        if(!inherits(object, "glm")){
          s2 <- sum(object$weights*object$residuals^2)/object$df.residual
          se.psi <- sqrt(diag(R)*s2)[idV]
        } else {
          s2<- if(object$fam$fam%in%c("poisson","binomial")) 1 else object$deviance/object$df.residual
          se.psi <- sqrt(diag(R)*s2)[idV]
        }
      object[["psi"]][,"St.Err"] <- se.psi
      }
    }
    if("segmented.default" == as.character(object$call)[1]){
      summ <- c(summary(object, ...), object["psi"])
      summ[c("it","epsilon")]<-object[c("it","epsilon")]
      #v<-try(vcov(object), silent=TRUE)
      #if(class(v)!="try-error") v<-sqrt(diag(v))
      return(summ)
      }
    if("lm"%in%class(object) && !"glm"%in%class(object)){
      summ <- c(summary.lm(object, ...), object["psi"])
      summ$Ttable<-summ$coefficients
      if(var.diff){
            #modifica gli SE
            Qr <- object$qr
            p <- object$rank #n.parametri stimati
            p1 <- 1L:p
            inv.XtX <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
            X <- qr.X(Qr,FALSE) 
            attr(X, "assign") <- NULL
            K<-length(unique(object$id.group)) #n.gruppi (=n.psi+1)
            dev.new<-tapply(object$residuals, object$id.group, function(.x){sum(.x^2)})
            summ$df.new<-tapply(object$residuals, object$id.group, function(.x){(length(.x)-eval(parse(text=p.df)))})
            if(any(summ$df.new<=0)) stop("nonpositive df when computig the group-specific variances.. reduce 'p.df'?", call. = FALSE)
            summ$sigma.new<-sqrt(dev.new/summ$df.new)
            sigma.i<-rowSums(model.matrix(~0+factor(object$id.group))%*%diag(summ$sigma.new))
            var.b<-inv.XtX%*%crossprod(X*sigma.i)%*%inv.XtX #sqrt(rowSums((X %*% V) * X))
            dimnames(var.b)<-dimnames(summ$cov.unscaled)
            summ$cov.var.diff<-var.b
            summ$Ttable[,2]<-sqrt(diag(var.b))
            summ$Ttable[,3]<-summ$Ttable[,1]/summ$Ttable[,2]
            summ$Ttable[,4]<- 2 * pt(abs(summ$Ttable[,3]),df=object$df.residual, lower.tail = FALSE)
            dimnames(summ$Ttable) <- list(names(object$coefficients)[Qr$pivot[p1]],
                  c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
            }
            if(!is.null(.vcov)){
              summ$Ttable[,2]<-sqrt(diag(.vcov))
              summ$Ttable[,3]<-summ$Ttable[,1]/summ$Ttable[,2]
              summ$Ttable[,4]<- 2 * pt(abs(summ$Ttable[,3]),df=object$df.residual, lower.tail = FALSE)
            #dimnames(summ$Ttable) <- list(names(object$coefficients)[Qr$pivot[p1]], c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
            }
      
      summ$Ttable[idU,4]<-NA
      summ$Ttable<-summ$Ttable[-idV,] 
      summ[c("it","epsilon","conv.warn")]<-object[c("it","epsilon","id.warn")]
      summ$n.boot<-length(na.omit(object$psi.history$all.ss))
      summ$var.diff<-var.diff
      summ$short<-short
      class(summ) <- c("summary.segmented", "summary.lm")
      return(summ)
      }
    #if("glm"%in%class(object)){
    if(inherits(object, "glm")){
      summ <- c(summary.glm(object, ...), object["psi"])
      summ$Ttable<-summ$coefficients[-idV,]
      summ$Ttable[idU,4]<-NA
      summ[c("it","epsilon","conv.warn")]<-object[c("it","epsilon","id.warn")]
      summ$n.boot<-length(na.omit(object$psi.history$all.ss))
      summ$short<-short
      class(summ) <- c("summary.segmented", "summary.glm")
      return(summ)}
    if("Arima"%in%class(object)){
      #da controllare
      coeff<-object$coef
      v<-sqrt(diag(object$var.coef))
      Ttable<-cbind(coeff[-idV],v[-idV],coeff[-idV]/v[-idV])
      colnames(Ttable)<-c("Estimate","Std. Error","t value")
      object$Ttable<-Ttable
      object$short<-short
      summ<-object
      summ[c("it","epsilon","conv.warn")]<-object[c("it","epsilon","id.warn")]
      summ$n.boot<-length(na.omit(object$psi.history$all.ss))
      class(summ) <- c("summary.segmented", "summary.Arima")
      return(summ)}
}

