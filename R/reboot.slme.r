reboot.slme <-function(fit, B=10, display=FALSE, metodo=1, frac=1, it.max=6, it.max.b=5, seed=NULL, start=NULL, msg=TRUE){
  #metodo: viene passato alla funzione logL. Se 1 la logL che viene calcolata e' quella della componente
  #   fit$lme.fit.noG, namely the logLik from the lme fit without the G variables..
  #bootRestart for slme4
  #fit: un oggetto di classe "segmented.lme" (anche proveniente da un altra "bootsegMix" call)
  #frac: size of the boot resample..
  #start : un vettor con i nomi (se non fornito gli starting values sono presi da fit)
  #-----------------------
  extract.psi<-function(obj){
    #questa funzione restituisce i "kappa", ovvero i coeff di psi..
    nomiG<-obj$namesGZ$nomiG
    b<-fixef(obj[[1]])[c("G0",nomiG)]
    b
  }
  
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
  #---------
  startKappa00<-extract.psi(fit)[1]
  Z <- fit$Z #segmented covariate
  rangeZ<-quantile(Z, c(.05,.95), names=FALSE)
  #quanti soggetti? Attenzione se ci sono nested re, sotto non funziona, o meglio da i livelli del outermost group
  
  #idLevels <- levels(fit$lme.fit$groups[,ncol(fit$lme.fit$groups)])
  #N<- length(idLevels)
  
  newData<-fit$lme.fit$data
  nomeRispo<-all.vars(formula(fit$lme.fit))[1]
  #AGGIUSTA la risposta
  newData[,nomeRispo]<-newData[,nomeRispo] + fit$Off
  
  nome.id <-names(fit$lme.fit$groups)[ncol(fit$lme.fit$groups)] #name of the innermost grouping variable 
  newData[, nome.id]<- factor(newData[, nome.id])
  var.id<-newData[, nome.id]
  idLevels<-levels(var.id)
  N<- length(idLevels)
  
  o.b<-fit$boot.call
  #old:    start.psi<-extract.psi(fit)
  #old:    est.psi<-start.psi["G0"]
  #old:    call.b<-update(object=fit, obj=o.b, data=newD, psi=est.psi, display=FALSE, evaluate=FALSE)
  call.b<-update(object=fit, obj=o.b, data=newD, it.max=it.max.b,
                 start=list(kappa0=startKappa0,kappa=startingKappa), display=FALSE, evaluate=FALSE)
  
  call.b$random <- fit$randomCALL
  
  o.ok<-update.lme.call(o.b, fixed.=paste(nomeRispo,"~."), evaluate=FALSE)
  #o.ok<-update.lme.call(o.b, fixed.=y~., evaluate=FALSE)
  #mycall$data=quote(gh)
  #o.ok<-update.lme.call(o.b, fixed.=y~.,evaluate=FALSE)
  #old:    call.ok<-update(object=fit, obj=o.ok, data=newData, psi=est.psi.b, display=FALSE, evaluate=FALSE)
  #o.ok$fixed<- update.formula(o.ok$fixed, paste(nomeRispo,"~."))
  
  call.ok<-update(object=fit, obj=o.ok, data=newData, it.max=it.max,
                  start=list(kappa0=startKappa0.b,kappa=startingKappa.b), display=FALSE, evaluate=FALSE)
  
  
  call.ok$n.boot <- call.b$n.boot<-0
  call.ok$control <- call.b$control<-quote(seg.control(display=FALSE))
  all.L<-all.psi<-NULL
  it<-0
  L0<-L.orig<-logLik(fit$lme.fit.noG)# logL(fit, metodo=metodo)
  if(display){
    flush.console()
    cat("original data:", 0, "  logLik =", formatC(as.numeric(L.orig), 3, format = "f"),"   psi parms:", formatC(extract.psi(fit),4,format="f"),"\n")
  }
  if(is.null(start)){
    startingKappa<-extract.psi(fit)
    startKappa0<- startingKappa[1]
    startingKappa<-startingKappa[-1]
    nomiKappa<-names(startingKappa)
    nomiKappa<-sapply(strsplit(nomiKappa, "G\\."),function(x)x[2])
    names(startingKappa) <- nomiKappa
  } else {
    nomiG<-sapply(strsplit(fit$namesGZ$nomiG, "G\\."),function(x)x[2])
    if(length(intersect(names(start), c("G0", nomiG)))!=length(start)) stop("'start' should include all the changepoint parameters")
    startKappa0<-start["G0"]
    startingKappa<-start[-which("G0"%in%names(start))]
    nomiKappa<-names(startingKappa)
  }
  if(is.null(seed)) seed<-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
  if(!is.numeric(seed)) stop(" 'seed' is not numeric")
  set.seed(seed)
  
  #browser()
  for(i in seq(B)){
    #build the boot sample
    #idx<-sample(N, replace=TRUE)
    #idx<-sample(1:N, size=trunc(N*frac), replace=TRUE)
    idx<-sample(idLevels, size=trunc(N*frac), replace=TRUE)
    
    newD <- do.call("rbind",lapply(idx, function(x)newData[newData[,nome.id]==x,]))
    newD$y.b<- newD[,nomeRispo]
    
    #       r<-list(newD=newD, call.b=call.b)
    #       return(r)
    
    #-->>       CAMBIA STARTING VALUE in call.b
    if(startKappa0>=rangeZ[2] | startKappa0<=rangeZ[1] ) startKappa0<- jitter(startKappa00,factor=5) #sum(rangeZ)/2
    
    fit.b<-try(suppressWarnings(eval(call.b)), silent=TRUE) #envir=newD) 
    if(!is.list(fit.b)){
      #        fit.b<-NULL
      it.b<-0
      while(!is.list(fit.b)){
        idx<-sample(idLevels, size=trunc(N*frac), replace=TRUE)
        newD <- do.call("rbind",lapply(idx, function(x)newData[newData[,nome.id]==x,]))
        newD$y.b<- newD[,nomeRispo]
        startKappa0<- jitter(startKappa00,factor=5)
        fit.b<-try(suppressWarnings(eval(call.b)), silent=TRUE) #envir=newD)
        it.b<-it.b+1
        if(it.b>=10) break
      }
    }
    if(is.list(fit.b)){
      #old: start.psi.b<-extract.psi(fit.b)
      #old: est.psi.b<-start.psi.b["G0"]
      startingKappa.b<-extract.psi(fit.b)
      startKappa0.b<- startingKappa.b[1]
      startingKappa.b<-startingKappa.b[-1]
      #NB "nomiKappa" dovrebbero essere sempre gli stessi
      names(startingKappa.b) <- nomiKappa
      fit.ok<-try(suppressWarnings(eval(call.ok)), silent=TRUE) # data=newData)
      #L1<-if(is.list(fit.ok)) logL(fit.ok, metodo=metodo) else (-Inf)
      #22/05/18 aggiunto un altro tentativo... ho notato che l'insuccesso pu? dipendere dagli starting value..
      if(!is.list(fit.ok)){
        call.ok$start<-NULL
        fit.ok<-try(suppressWarnings(eval(call.ok)), silent=TRUE)
      }
      L1<-if(is.list(fit.ok)) as.numeric(logLik(fit.ok)) else (-Inf)
    } else {
      stop("the bootstrap fit is unsuccessful")
    }
    if(L0<L1) {
      fit<-fit.ok
      L0<-L1
    }
    all.psi[length(all.psi)+1]<-est.psi<-extract.psi(fit)["G0"]
    all.L[length(all.L)+1]<-L.ok<-max(L0,L1)
    it<-it+1
    if(display){
      flush.console()
      ll<-if(it<10) "  logLik =" else " logLik ="
      cat("boot resample:", it, ll, formatC(L.ok, 3, format = "f"),"   psi parms:", formatC(extract.psi(fit),4,format="f"),"\n")
    }
    startingKappa<-extract.psi(fit)
    startKappa0<- startingKappa[1]
    startingKappa<-startingKappa[-1]
    nomiKappa<-names(startingKappa)
    nomiKappa<-sapply(strsplit(nomiKappa, "G\\."),function(x)x[2])
    names(startingKappa) <- nomiKappa
  } #end boot replicates
  fit$history.boot.restart<-cbind(b=1:length(all.psi),psi=all.psi, logL=all.L)
  fit$seed<-seed
  #r<-list(seg.lme.fit=fit, history=cbind(b=1:length(all.psi),psi=all.psi, logL=all.L) )
  if(msg) cat(" New solution(s) found:", length(unique(all.psi))-1, "\n")
  fit
}


