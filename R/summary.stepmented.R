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
    nomiPsi<- sub("V", "psi", nomiV)
    idU <-match(nomiU,names(coef(object)[!is.na(coef(object))]))
    idV <-match(nomiPsi,names(coef(object)[!is.na(coef(object))]))
    #beta.c<- coef(object)[nomiU]
    #per metodo default.. ma serve????
    if("stepmented.default" == as.character(object$call)[1]){
      summ <- c(summary(object, ...), object["psi"])
      summ[c("it","epsilon")]<-object[c("it","epsilon")]
      #v<-try(vcov(object), silent=TRUE)
      #if(class(v)!="try-error") v<-sqrt(diag(v))
      return(summ)
    }
    #browser()
    
    VAR <- if(!is.null(.vcov)) .vcov else vcov(object,...)
    se <- sqrt(diag(VAR))
    object$psi[,"St.Err"] <- se[nomiPsi] 
    
    #if("lm"%in%class(object) && !"glm"%in%class(object)){
    if(inherits(object, "lm") && !inherits(object, "glm")){
      #object$rank include i psi, mentre object$qr$rank no.
      #Affinche' summary.lm() funzioni, e' necessario che object$rank non tenga conto del numero di psi..
      #quindi qua (e anche nei lm sopra) modifichiamo il valore di rank..
      #NB: questo problema NON si presenta se sono state usate le funzioni stepmented.* in cui anche object$qr$rank tiene gia' conto
      # dei psi (perche' hanno stimato il modello con le variabili W per cercare di ottenere una qualche misura del se)
      object$rank <- object$qr$rank #object$rank - nrow(object$psi)
      # summ <- c(suppressWarnings(summary.lm(object, ...)), object["psi"])
      # summ$Ttable <-summ$coefficients
      # b <- coef(object, FALSE)
      # b<-b[b!=0]
      # summ <- list(Ttable=matrix(NA, length(b), 4, dimnames = list(names(b),c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))),
      #              psi=object[["psi"]], sigma=sigma(object), call=object$call, df=c(length(coef(object)), object$df.residual, length(coef(object)) ) )
      # summ$Ttable[,"Estimate"] <- b 
      # summ$Ttable[,"Std. Error"] <- se[1:length(b)] #se[rownames(summ$coefficients)]
      # summ$Ttable[,"t value"] <- summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"] 
      # summ$Ttable[,"Pr(>|t|)"] <- 2*pt(abs(summ$Ttable[,"t value"]), df=object$df.residual, lower.tail = FALSE) # summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"]
      
      summ <- c(suppressWarnings(summary.lm(object, ...)), object["psi"])
      summ$Ttable <- summ$coefficients
      summ$Ttable[, "Std. Error"] <- se[rownames(summ$coefficients)]
      summ$Ttable[, "t value"] <- summ$Ttable[, "Estimate"]/summ$Ttable[,"Std. Error"]
      summ$Ttable[, "Pr(>|t|)"] <- 2 * pt(abs(summ$Ttable[,"t value"]), df = object$df.residual, lower.tail = FALSE)
      
      
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
      summ$Ttable[idU,4]<-NA
      if(all(!is.na(idV))) summ$Ttable<-summ$Ttable[-idV,] 
      summ[c("it","epsilon","conv.warn")]<-object[c("it","epsilon","id.warn")]
      summ$n.boot<-length(na.omit(object$psi.history$all.ss))
      summ$var.diff<-var.diff
      summ$short<-short
      summ$psi.rounded <- object$psi.rounded
      class(summ) <- c("summary.stepmented", "summary.lm")
      return(summ)
    }
    if(inherits(object, "glm")){
      #browser()
      #23/4/24 mi sono reso conto che con gaussian GLM viene stampato "t-value" e non z-value... 
      #     Per cui piuttostoche i nomi, metto gli indici delle colonne..
      object$rank <- object$qr$rank
      summ <- c(suppressWarnings(summary.glm(object, ...)), object["psi"])
      summ$Ttable <-summ$coefficients
      summ$Ttable[,"Std. Error"] <- se[rownames(summ$coefficients)]
      # summ$Ttable[,"z value"] <- summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"] 
      # summ$Ttable[,"Pr(>|z|)"] <- 2*pnorm(abs(summ$Ttable[,"z value"]), lower.tail = FALSE) # summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"]
      summ$Ttable[,3] <- summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"] 
      summ$Ttable[,4] <- if(object$family$family=="gaussian") 2*pt(abs(summ$Ttable[,3]), df=object$df.residual, lower.tail = FALSE) else 2*pnorm(abs(summ$Ttable[,3]), lower.tail = FALSE) # summ$Ttable[,"Estimate"]/summ$Ttable[,"Std. Error"]
      
      summ$Ttable[idU,4]<-NA
      if(all(!is.na(idV))) summ$Ttable<-summ$Ttable[-idV,] 
      summ[c("it","epsilon","conv.warn")]<-object[c("it","epsilon","id.warn")]
      summ$n.boot<-length(na.omit(object$psi.history$all.ss))
      summ$short<-short
      summ$psi.rounded <- object$psi.rounded
      class(summ) <- c("summary.stepmented", "summary.glm")
      return(summ)
    }
    if("Arima"%in%class(object)){
      stop("stepmented arima model not allowed")
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

