`summary.stepmented` <-
function(object, short=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, ...){
    
    if(!is.null(.vcov)) var.diff<-FALSE
    if(var.diff && length(object$nameUV$Z)>1) {
      var.diff<-FALSE
      warning(" 'var.diff' set to FALSE with multiple segmented variables", call.=FALSE)
      }
    
    #browser()
    
    nomiU<-object$nameUV$U
    nomiV<-object$nameUV$V
    nomiPsi<- gsub("V", "psi", nomiV)
    idU <-match(nomiU,names(coef(object)[!is.na(coef(object))]))
    idV <-match(nomiPsi,names(coef(object)[!is.na(coef(object))]))
    beta.c<- coef(object)[nomiU]
    #per metodo default.. ma serve????
    if("stepmented.default" == as.character(object$call)[1]){
      summ <- c(summary(object, ...), object["psi"])
      summ[c("it","epsilon")]<-object[c("it","epsilon")]
      #v<-try(vcov(object), silent=TRUE)
      #if(class(v)!="try-error") v<-sqrt(diag(v))
      return(summ)
      }
    if("lm"%in%class(object) && !"glm"%in%class(object)){
      summ <- c(summary.lm(object, ...), object["psi"])
      summ$Ttable <-summ$coefficients
      
      summ$Ttable[,"Std. Error"] <- sqrt(diag(vcov(object)))
      summ$Ttable[,"t value"] <- summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"] 
      summ$Ttable[,"Pr(>|t|)"] <- 2*pt(abs(summ$Ttable[,"t value"]), df=object$df.residual, lower.tail = FALSE) # summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"]
      
      #browser()
      
      if(var.diff){
        stop("not allowed")
        # modifica gli SE
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
      }
          summ$Ttable[idU,4]<-NA
          summ$Ttable<-summ$Ttable[-idV,] 
          summ[c("it","epsilon","conv.warn")]<-object[c("it","epsilon","id.warn")]
          summ$n.boot<-length(na.omit(object$psi.history$all.ss))

          summ$var.diff<-var.diff
          summ$short<-short
          summ$psi.rounded <- object$psi.rounded
          class(summ) <- c("summary.stepmented", "summary.lm")
          return(summ)
    }
    if(inherits(object, "glm")){
      summ <- c(summary.glm(object, ...), object["psi"])
      summ$Ttable<-summ$coefficients[-idV,]
      #cat("HAi modiifcato il calcolo degli SE comesta fatto sopra per 'lm'? \n")
      summ$Ttable[idU,4]<-NA
      summ[c("it","epsilon","conv.warn")]<-object[c("it","epsilon","id.warn")]
      summ$n.boot<-length(na.omit(object$psi.history$all.ss))
      summ$short<-short
      summ$psi.rounded <- object$psi.rounded
      class(summ) <- c("summary.stepmented", "summary.glm")
      return(summ)
      }
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
      summ$psi.rounded <- object$psi.rounded
      class(summ) <- c("summary.stepmented", "summary.Arima")
      return(summ)
    }
}

