vcov.segmented.lme <-function(object, B=0, ret.b=FALSE, ...){
  bootNP<-function(fit, B=50, seed=NULL, it.max.b=6){
    #Non parametric boot for slme4
    #fit: un oggetto di classe "segmented.lme"
    #-----------------------
    update.lme.call<-function (old.call, fixed., ..., evaluate=FALSE) {
      call <- old.call
      extras <- match.call(expand.dots = FALSE)$...
      if (!missing(fixed.)) call$fixed <- update.formula(call$fixed, fixed.)
      if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
          call <- c(as.list(call), extras[!existing])
          call <- as.call(call)
        }
      }
      if (evaluate) eval(call, parent.frame()) else call
    }
    #---------
    if(is.null(B) || B<=0) stop("'B>0' is requested")
    newData<-fit$lme.fit$data
    rnfGrp<- fit$lme.fit.noG$groups
    if(ncol(rnfGrp)>1) warning("the innermost grouping variable is used", call. = FALSE)
    nome.id <-names(rnfGrp)[ncol(rnfGrp)] #name of the innermost grouping variable    
    var.id<-newData[, nome.id]
    #idLevels<-levels(var.id)
    idLevels<-levels(fit$lme.fit$groups[[ncol(rnfGrp)]])
    N<-nlevels(fit$lme.fit$groups[[ncol(rnfGrp)]]) #n. of "subjects"
    
    nomeRispo<-all.vars(formula(fit$lme.fit))[1]
    #AGGIUSTA la risposta
    newData[,nomeRispo]<-newData[,nomeRispo] + fit$Off
    o.b<-fit$boot.call
    call.b<-update(object=fit, obj=o.b, data=newD, it.max=it.max.b,
                   start=list(kappa0=startKappa0,kappa=startingKappa), display=FALSE, evaluate=FALSE)
    call.b$random <- fit$randomCALL
    startingKappa<-extract.psi(fit)
    startKappa0<- startingKappa[1]
    startingKappa<-startingKappa[-1]
    nomiKappa<-names(startingKappa)
    nomiKappa<-sapply(strsplit(nomiKappa, "G\\."),function(x)x[2])
    names(startingKappa) <- nomiKappa
    
    est<-fixef(fit$lme.fit)
    se<-sqrt(diag(vcov(fit$lme.fit)))
    fitt<-fitted.segmented.lme(fit, level=0)
    COEF<-SE<-matrix(,B,length(est))
    FIT<-matrix(,B,length(fitt))
    if(!is.null(seed)) set.seed(seed)
    
    for(i in seq(B)){
      #build the boot sample
      #idx<-sample(N, replace=TRUE)
      #idx<-sample(1:N, size=N, replace=TRUE)
      #idx<-levels(fit$lme.fit$groups[[1]])[idx]
      #newD<-do.call("rbind",lapply(idx, function(x)newData[newData$id==x,]))
      #newD$y.b<- newD$y
      idx<-sample(idLevels, size=N, replace=TRUE) 
      newD <- do.call("rbind",lapply(idx, function(x)newData[newData[,nome.id]==x,]))
      newD$y.b<- newD[,nomeRispo]
      
      fit.b<-try(suppressWarnings(eval(call.b)), silent=TRUE) #envir=newD)
      if(is.list(fit.b)){
        Tt<- summary(fit.b[[1]])$tTable #usa summary.lme
        COEF[i,]<-Tt[,1] #coef
        SE[i,]<-Tt[,2] #se
        FIT[i,]<-fitted.segmented.lme(fit.b, level=0)
      }
    }
    r<-list(coef=rbind(est,COEF),se=rbind(se,SE), fitted=rbind(fitt,FIT))
    r
  }
  #--------------------
  extract.psi<-function(obj){
    #questa funzione restituisce i "kappa", ovvero i coeff di psi..
    nomiG<-obj$namesGZ$nomiG
    b<-fixef(obj[[1]])[c("G0",nomiG)]
    b
  }
  #--------------------
  opz<-list(...)
  opz$fit <- object
  opz$B <- B
  if(B<=0) {
    r <- object$lme.fit$varFix
  } else {
    #browser()
    r <- do.call(bootNP, opz) #bootNP(object, B= B, seed=opz$seed, opz$it.max.b=6)
    if(!ret.b) r<- var(r$coef[-1,])
  }
  return(r)
}
