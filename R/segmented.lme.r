segmented.lme <- function(obj, seg.Z, psi, npsi=1, fixed.psi=NULL, control = seg.control(), model = TRUE,
                          z.psi=~1, x.diff=~1, 
                          random=NULL, #una lista quale 'list(id=pdDiag(~1+x+U+G0))'
                          random.noG=NULL, #una lista senza G0. Se NULL viene aggiornata la formula di random escludendo "G0"
                          start.pd=NULL, #una matrice come starting value
                          psi.link=c("identity","logit"), 
                          #nq=0, 
                          #adjust=0,
                          start=NULL, #*named* list list(delta0, delta, kappa) and the 'delta' component, dovrebbe essere anche
                          #nominata con i nomi delle variabili in x.diff
                          data,
                          fixed.parms=NULL,...){ #a *named* vector meaning the coefficients to be mantained fixed during the estimation
                          #, tol=0.0001, it.max=10, display=FALSE){
  #control = list(niterEM = 0, optimMethod = "L-BFGS-B")
  #method = "ML"
  ################################################################################
  
  #require(nlme)
  adj.psi <- function(psii, LIM) {
    pmin(pmax(LIM[1, ], psii), LIM[2, ])
  }
  
  newData<-aa<-betaa<-fn1<-kappa1<-NULL
  tol <- control$toll
  it.max <- control$it.max
  display <- control$visual
  n.boot <- control$n.boot
  alpha <- control$alpha
  if(is.null(alpha)) alpha<- max(.05, 1/obj$dims$N)
  if(length(alpha)==1) alpha<-c(alpha, 1-alpha)
  
  adjust=0 #ho rimosso dagli argomenti adjust=0, pero' devo ancora vederlo bene..
  
  
  psi.link<-match.arg(psi.link)
  logit<-function(xx,a,b){log((xx-a)/(b-xx))}
  inv.logit<-function(xx,a,b){((a+b*exp(xx))/(1+exp(xx)))}
  
  #obj is the lme fit or simply its call
  #random: a list with a formula for the cluster variable 'id' and standard linear variables and "U" and "G0" meaning
  #     random effects for the difference in slope and changepoint parameters. If it.max=0 the breakpoint is not estimated and
  #     the formula should not include the term "G0".
  #random = list(id=pdBlocked(list(pdDiag(~1+x), pdSymm(~U+G0-1))))
  #random = list(id=pdBlocked(list(pdSymm(~1+x), pdSymm(~U+G0-1))))
  #random=list(id=pdDiag(~1+weeks+U+G0))
  #random=list(id=pdSymm(~1+weeks+U+G0))
  #
  #Problemi: se control?
  #control = list(msVerbose = FALSE, niterEM = 100, opt = "optim")
  #
  #nq: no. obs che consentono di "invalidare" la stima del breakpoints.
  # Ovvero se nq=0, gli \hat{\psi}_i sono annullati se \hat{\psi}_i<=min(Z_i) o \hat{\psi}>=max(z_i)
  #        se nq>0 gli \hat{\psi}_i sono annullati se \hat{\psi}_i<=min(sort(z)[1:nq]) o \hat{\psi}>= max(rev(z)[1:nq]
  #adjust valore numerico (0,1,2).
  #   Se 0 i psi_i vengono stimati "normalmente" e alla convergenza al vettore numerico dei psi viene assegnato un
  #   vettore di attributi che serve ad etichettare se il breakpoint ? plausibile o meno (secondo il valore di nq)
  #   Se 1 i psi ottenuti alla fine dell'algoritm vengono aggiustati secondo il valore di nq. Ad es., se nq=1 il breakpoint
  #   immediatamente prima del max (o dopo il min) vengono forzati al min/max e cos? sono di fatto annullati; naturalmente il
  #   modello ? ristimato secondo  i nuovi psi. Se 2 l'aggiustamento viene fatto durante l'algoritmo..
  #---------------------
  reboot.slme <-function(fit, B=10, display=FALSE, break.boot=B, metodo=1, frac=1, it.max=6, it.max.b=5, seed=NULL, start=NULL, msg=TRUE){
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
                    start=list(kappa0=startKappa0.b, kappa=startingKappa.b), display=FALSE, evaluate=FALSE)
    
    
    #call.ok$n.boot <- call.b$n.boot<-0
    call.ok$control <- call.b$control<-quote(seg.control(display=FALSE, n.boot=0))
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
    #if(is.null(seed)) seed<-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
    if(is.null(seed)){
      mY <- mean(newData[,nomeRispo])
      sepDec<-if(options()$OutDec==".") "\\." else "\\,"
      vv <- strsplit(paste(strsplit(paste(mY), sepDec)[[1]], collapse=""),"")[[1]]
      vv<-vv[vv!="0"]
      vv=na.omit(vv[1:5])
      seed <-eval(parse(text=paste(vv, collapse="")))
      set.seed(seed)
    } else {
      if(is.na(seed)) {
        seed <-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
        set.seed(seed)
      } else {
        if(!is.numeric(seed)) stop(" 'seed' is not numeric") else set.seed(seed)
      }
    }  
    

    #browser()
    n.boot.rev<- 3
    alpha1<-alpha[1]
    for(i in seq(B)){
      diff.selected.ss <- rev(diff(na.omit(all.L)))
      if(length(diff.selected.ss)>=(n.boot.rev-1) && all(round(diff.selected.ss[1:(n.boot.rev-1)],6)==0)){
        #qpsi<-sapply(1:ncol(Z),function(i)mean(est.psi0[i]>=Z[,i]))
        qpsi<- mean(startKappa0>Z)
        qpsi<-ifelse(abs(qpsi-.5)<.1, alpha1, qpsi)
        alpha1<-1-alpha1
        #est.psi0<-sapply(1:ncol(Z),function(i)quantile(Z[,i],probs=1-qpsi[i],names=FALSE))
        startKappa0 <- quantile(Z, probs=1-qpsi, names=FALSE)
      }
      
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
      
      #conta i valori ss uguali per fermarsi prima..
       asss<-na.omit(all.L)
       if(length(asss)>break.boot){
         if(all(rev(round(diff(asss),6))[1:(break.boot-1)]==0)) break
       }
      
    } #end boot replicates
    #============================================================================================
    fit$history.boot.restart<-cbind(b=1:length(all.psi),psi=all.psi, logL=all.L)
    fit$seed<-seed
    #r<-list(seg.lme.fit=fit, history=cbind(b=1:length(all.psi),psi=all.psi, logL=all.L) )
    if(msg) cat(" New solution(s) found:", length(unique(all.psi)), "\n")
    fit
  }
  #------------------
  fn.re<-function(obj){
    #restituisce un array n x n.ranef x terms
    #   n e' il n. totale delle misurazioni..
    #   n.ranef e' il n. dei random effects (tipicamente e' 1, >1 con nested..)
    #   terms e' il n. dei termini coinvolti nei random effects (ad es., intercept, x ..)
    ro<-ranef(obj)
    n.levels<- ncol(obj$groups) #n. dei livelli casuali (ad es., se nested..)
    if(n.levels<=1) {
      ro<-list(ro)
      names(ro)<-names(obj$groups)
    }
    nomi.levels<-names(obj$groups) #nomi degli effetti casuali names(ranef(obj))
    n.terms<-sapply(ro, ncol)
    nomiTermini<- unique(as.vector(unlist(sapply(ro, colnames))))
    tutti<-array(0, c(nrow(obj$groups), ncol(obj$groups), max(n.terms)), dimnames=list(NULL, names(obj$groups), nomiTermini))
    for(nome in nomiTermini){
      for(j in nomi.levels){
        if(nome %in% names(ro[[j]])){
          for(i in unique(obj$groups[,j])) tutti[obj$groups[,j]==i,j,nome] <- ro[[j]][rownames(ro[[j]])==i, nome]
        }
      }
    }
    tutti
  }
  #------------------
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
  #---------------------------------------------------------------------------
  f.pd<-function(obj){
    #dato un modello lme 'obj' restituisce una matrice pdMat che deve essere utilizzata come componente random
    #   nelle call "call.ok$random<-list(id=pd)"
    pdClasse<-class(obj$modelStruct$reStruct[[1]])[1]
    if(pdClasse=="pdBlocked"){ #assumiamo solo 2 blocchi..(? un LIMITE, ma ? facile generalizzare..)
      start.v<-unlist(lapply(obj$modelStruct$reStruct[[1]], function(z){as.numeric(z)}))
      cl1<-class(obj$modelStruct$reStruct[[1]][[1]])[1]
      cl2<-class(obj$modelStruct$reStruct[[1]][[2]])[1]
      fo1<-attr(obj$modelStruct$reStruct[[1]][[1]],"formula")
      fo2<-attr(obj$modelStruct$reStruct[[1]][[2]],"formula")
      no1<-attr(obj$modelStruct$reStruct[[1]][[1]],"Dimnames")[[1]]
      no2<-attr(obj$modelStruct$reStruct[[1]][[2]],"Dimnames")[[1]]
      pd<-pdBlocked(start.v, pdClass = c(cl1,cl2), nam = list(no1, no2), form=list(fo1, fo2))
    } else {
      fo<-attr(obj$modelStruct$reStruct[[1]],"formula")
      pd <- pdMat(as.numeric(obj$modelStruct$reStruct[[1]]), form = fo, pdClass = pdClasse)
    }
    pd}
  #---------------------------------------------------------------------------
  ###
  #browser()
  h <- control$h 
  if(!(is.call(obj) || class(obj)[1]=="lme")) stop(" 'obj' should be a lme fit or a lme call")
  if(missing(psi) && it.max==0) stop("Please supply 'psi' with 'it.max=0'")
  
  if(is.call(obj)) {
    my.call  <- obj
    datacall <- deparse(obj$data)
    if(is.null(random)) random<-eval(obj$random)      
  } else {
    my.call <- obj$call
    datacall<- deparse(obj$call$data)
    if(is.null(random)) random<-eval(obj$call$random)  
  }
  #my.call<-if(is.call(obj)) obj else obj$call
  #datacall<- if(is.call(obj)) deparse(obj$data) else deparse(obj$call$data)
  #if(is.null(random)) {random<- if(is.call(obj)) eval(obj$random) else eval(obj$call$random) }
  randomCALL<-random
  G0random<- sapply(random, function(.x) "G0" %in% all.vars(attr(.x, "formula")))
  if(it.max==0 && !any(G0random)) stop("'G0' in the random part is meaningless with 'it.max=0'")
  #    name.group<-nameRandom<-names(random)
  
  #    if(is.null(random)) {
  #      # A CHE SERVE????????????????
  #      random=list(
  #          id=pdMat(as.numeric(obj$modelStruct$reStruct[[1]]),
  #          form=attr(obj$modelStruct[[1]][[1]],"formula"),
  #          pdClass=class(obj$modelStruct$reStruct[[1]])[1]))
  
  #     randomCALL<- if(is.call(obj)) obj$random else obj$call$random
  #     } else {
  #	      randomCALL<- random
  #    }
  
  if (!is.null(random)) {
    if (is.list(random)) {
      nameRandom <- names(random) #nomi dei fattori id
      if(is.null(nameRandom)) stop("random argument must be a *named* list.")
      else if(sum(nameRandom == "")) stop("all elements of random list must be named")
    } else stop("random effects should be specified as named lists")
    random.vars <- c(unlist(lapply(random, function(x) all.vars(formula(x)))), nameRandom)
    names(random.vars)<-NULL #per evitare casini.. spesso i nomi erano le variabili stesse..
  } else random.vars <- NULL
  
  J<-length(random)
  
  #if(missing(Z) && missing(seg.Z)) stop(" 'Z' or 'seg.Z' should be provided")
  #name.Z<-if(missing(seg.Z)) deparse(substitute(Z)) else all.vars(seg.Z)
  
  if(missing(seg.Z)) stop(" 'seg.Z' should be provided")
  name.Z<- all.vars(seg.Z)
  if(length(name.Z)>1) stop("segmented.lme works with 1 breakpoint only")
  
  allNOMI<-unique(c(name.Z, all.vars(my.call$fixed), random.vars, all.vars(z.psi), all.vars(x.diff)))
  formTUTTI<-as.formula(paste("~.+", paste(allNOMI,collapse="+")))
  formTUTTI<-update.formula(my.call$fixed, as.formula(paste("~.+", paste(allNOMI,collapse="+"))))
  #U and G0 have not yet been defined
  formTUTTI<-update.formula(formTUTTI, .~.-U-G0)
  
  anyFixedG<-FALSE
  if(!is.null(fixed.parms)){
    name.fixed.butG0<-setdiff(names(fixed.parms),"G0") #nomi dei termini fissi escluso G0
    anyFixedG<-if(length(name.fixed.butG0)>=1) TRUE else FALSE #ci sono fixed coef nel submodel of psi?
    if(anyFixedG){
      formTUTTI<-update.formula(formTUTTI, as.formula(paste("~.+", paste(name.fixed.butG0,collapse="+"))))
    }
  }
  
  
  if(is.null(my.call$data)) stop("`obj' should include the argument `data'")
  if(missing(data)) {
    mf<-model.frame(formTUTTI, data=eval(my.call$data), na.action=na.omit)
  } else {
    mf<-model.frame(formTUTTI, data=data, na.action=na.omit)
  }
  
  
  
  #    if (length(allvars)) {
  #        mf$formula <- as.formula(paste(paste(deparse(gp$fake.formula, 
  #            backtick = TRUE), collapse = ""), "+", paste(allvars, 
  #            collapse = "+")))
  #        mf <- eval(mf, parent.frame())
  #    }
  
  #adesso si deve ordinare il dataframe..
  mf<-mf[order(mf[[nameRandom[J]]]),] 
  
  nomeRispo<-names(mf)[1]
  Rispo<-model.response(mf)
  #
  
  #browser()
  
  
  Z <- mf[[name.Z]]
  
  #limZ <- apply(Z, 2, quantile, names = FALSE, probs = c(alpha[1], alpha[2]))

  limZ <- as.matrix(quantile(Z, names = FALSE, probs = c(alpha[1], alpha[2])))
  
  min.Z<- min(limZ[,1])
  max.Z<- max(limZ[,1])
  
  
  
  
  
  
  
  #browser()
  
  if(!missing(psi)) {
    if(length(psi)>1) stop("segmented.lme works with 1 breakpoint only")
    if(psi<=min(limZ) || psi>=max(limZ)) stop("the provided psi is outside the range, see 'alpha' in seg.control()", call.=FALSE)
    }
  
  id <- mf[[nameRandom[J]]] #the innermost factor
  if(is.factor(id)) id <-factor(id, levels = unique(id)) 
  
  ni<- tapply(id, id, length) #vector of cluster sizes
  N<-length(ni)#n. of clusters (subjects)
  n<-length(id) #n. of total measurements
  
  id.x.diff<- FALSE
  id.z.psi <- FALSE
  #M.z.psi <- mf[all.vars(z.psi)] #
  #M.x.diff <- mf[all.vars(x.diff)] #
  
  M.z.psi <- model.matrix(z.psi, data = mf)
  if("(Intercept)"%in%colnames(M.z.psi)) M.z.psi<-M.z.psi[,-match("(Intercept)", colnames(M.z.psi)),drop=FALSE]
  M.x.diff <- model.matrix(x.diff, data = mf)
  if("(Intercept)"%in%colnames(M.x.diff)) M.x.diff<-M.x.diff[,-match("(Intercept)", colnames(M.x.diff)),drop=FALSE]
  
  fixed<-"U+G0" #fixed<-"U"
  nomiG<-NULL #se non ci sono explicative nel changepoint (se ci sono poi viene sovrascritto)
  namesGZ<-list(nameZ=name.Z)
  
  Offs.kappa<-0
  if(NCOL(M.z.psi)>0){
    id.z.psi <- TRUE
    Z.psi  <- data.matrix(M.z.psi)
    if(anyFixedG){
      if(!all(name.fixed.butG0 %in% colnames(M.z.psi))) stop("variable(s) in 'fixed.parms' should be included in 'z.psi'")
      Offs.kappa<-Fixed.z.psi<-drop(Z.psi[, name.fixed.butG0, drop=FALSE]%*% fixed.parms[name.fixed.butG0])
      Z.psi<-Z.psi[,setdiff(colnames(Z.psi), name.fixed.butG0), drop=FALSE]
    }
  if(ncol(Z.psi)>0){
      nomiG<-paste("G.",colnames(Z.psi),sep="") 
      namesGZ$nomiG<-nomiG
      fixed<-paste(fixed,paste(nomiG,collapse="+"),sep="+")
    } else {
      id.z.psi <- FALSE
    }
  } else { #se NCOL(M.z.psi)<=0
    if(anyFixedG) stop("variable(s) in 'fixed.parms' should be included in 'z.psi' ")
  }
  if(NCOL(M.x.diff)>0) {
    X.diff <- data.matrix(M.x.diff) 
    id.x.diff <- TRUE
    nomiUx<-paste("U.",colnames(M.x.diff),sep="")
    namesGZ$nomiUx<-nomiUx
    fixed<-paste(fixed,paste(nomiUx,collapse="+"),sep="+")
  }
  
  #==================================================================
  #Queste funzioni min1() e max1() restituiscono il "quasi" min o max
  # if(nq>0){
  #   min1<-function(x,na.rm=FALSE){x<-sort(x)[-(1:nq)];min(x,na.rm=na.rm)}
  #   max1<-function(x,na.rm=FALSE){x<-rev(x)[-(1:nq)];max(x,na.rm=na.rm)}
  # } else {
  #   min1<-min
  #   max1<-max
  # }
  # adjust<-max(min(adjust,2),0)  #solo 0,1,2 sono consentiti..
  # 
  # #==================================================================
  # 
  # min.Z<-min1(Z)
  # max.Z<-max1(Z)
  
  
  mf["U"]<- 1 #rep(1, n)
  #if(!is.null(obj$data)) my.dd<-cbind(obj$data,my.dd)
  #browser()
  #Qua ci possono essere 2 variabili di effetti casuali. Attenzione all'ordine.. il secondo!
  #if(name.group!="id") mf['id']<-mf[name.group] #costruisci un'altra variabile di clustering con il nome id
  #correzione per nested r.e: poich? id ? quello "giusto" (costruito prima), allora
  #
  mf['id']<-id #E' necessario costruire una nuova id con nome esattamente 'id'??!??!
  mf[name.Z]<- Z
  
  est.kappa0<-TRUE
  if("G0" %in% names(fixed.parms)) {
    est.kappa0<-FALSE
    kappa0<-kappa0Fixed<-fixed.parms["G0"]
  }

  if(est.kappa0){
    if(!is.null(start$kappa0)) {
      psi<-if(psi.link=="logit") inv.logit(start$kappa0,min.Z,max.Z) else start$kappa0
    }
    
    if(missing(psi)){
      #        formulaFix.Poly<-update.formula(my.call$fixed, paste("~.+",name.Z,"+",paste("I(",name.Z,"^2)",sep="")))
      #        obj2<-update.lme.call(my.call, fixed = formulaFix.Poly, data=mf, evaluate=TRUE)
      #        psi<- -fixed.effects(obj2)[name.Z]/(2*fixed.effects(obj2)[paste("I(",name.Z,"^2)",sep="")])
      psi<-tapply(Z, id, function(.x) sum(range(.x))/2)   
      if(any(psi <= min(Z))||any(psi>=max(Z))) stop("psi estimated by midvalues is outside the range") #the quadratic fit 
    }
  } else { #se e' fissato e quindi non devi stimarlo
    psi<- kappa0
  }
  
  
  #browser()
  
  
  psi.new <- psi #stime iniziali
  if(length(psi)!=1 && length(psi)!=N) stop("length(psi) has to be equal to 1 or n. of clusters")
  if(length(psi) == 1) {
    psi.new <- rep(psi.new, N) #subj-specific changepoints
  }
  psi.ex<-rep(psi.new, ni ) #length = N (n. tot obs)
  
  #----------------------------------------
  mf$U<- (Z-psi.ex)*(Z>psi.ex) #pmax(0, Z-psi.ex)
  formulaFix.noG<-update.formula(my.call$fixed, paste("~.+","U"))
  if(id.x.diff){
    Ux<- as.matrix(mf$U*X.diff)
    colnames(Ux)<-nomiUx
    mf<-cbind(mf,Ux) #$Ux<- my.dd$U*X.diff
    formulaFix.noG<-update.formula(my.call$fixed, paste(".~.+U+",paste(nomiUx,collapse="+"),sep=""))
  }
  #se vuoi assumere i psi fissi (it.max=0)
  if(it.max==0) {
    #aggiorna i random effects. Attenzione in tal caso random deve essere "U" ( o "1").
    #Se fosse "U+G0" darebbe errore perch? G0 non esiste
    #Oppure dovresti modificare la formula di random,
    #attr(random[[1]], "formula")<-update.formula(attr(random[[1]], "formula"), ~.-G0)
    formulaRand<-formulaRandOrig<-my.call$random
    call.ok<-update.lme.call(my.call, fixed = formulaFix.noG, random=random, data=mf, evaluate=FALSE)
    o<-eval(call.ok)
    return(o)
  } #end if(it.max=0)
  #---------------------------------------------------------------------------
  #should we fit a preliminary model? extract starting values
  start.delta0<-start$delta0
  if(id.x.diff) start.delta<-start$delta
  need.prelim<- (is.null(start.delta0) || (id.x.diff && is.null(start.delta)))
  
  if(need.prelim){
    random.noG <- random
    for(j in 1:J) attr(random.noG[[j]],"formula")<-update.formula(formula(random[[j]]), ~.-G0)
    o<-update.lme.call(my.call, fixed=formulaFix.noG, random=random.noG, data=mf, evaluate=TRUE)
    #o<-update.lme.call(my.call, fixed=formulaFix.noG, data=mf, evaluate=TRUE)
    delta0i<-unlist(coef(o)["U"]) #length= N
    if(id.x.diff) delta<-fixed.effects(o)[nomiUx] #length= n.1
  } else {
    delta0i<-if(length(start.delta0)==N) start.delta0 else rep(start.delta0,N)
    if(id.x.diff) delta<-start.delta[nomiUx]
  }
  
  start.kappa<-start$kappa
  
  eta.psi<-0
  
  if(id.z.psi) {
    if(is.null(start.kappa)) {
      kappa<- rep(0, ncol(Z.psi))
      names(kappa)<-nomiG
      eta.psi<-rep(0,nrow(Z.psi))
    } else {
      kappa<-start.kappa
      names(kappa)<-paste("G.",names(kappa),sep="")
      if((length(kappa)!=NCOL(M.z.psi)) || any(is.na(match(names(kappa), nomiG)))) stop("error in the names/length of start.kappa")
      eta.psi <- drop(Z.psi%*%kappa)
    }
  }
  #################################
  if(anyFixedG) eta.psi<- eta.psi + Offs.kappa
  #Offs.kappa<-data.matrix(mf[name.fixed.butG0])%*%fixed.parms[name.fixed.butG0]
  
  #-----------------------------------------------------------
  formulaFix<-update.formula(my.call$fixed, paste(".~.+",fixed))
  
  if(!est.kappa0) formulaFix<-update.formula(formulaFix, .~.-G0)
  formulaRand<-formulaRandOrig<-my.call$random
  minMax <- cbind(tapply(Z,id,min),tapply(Z,id,max)) #matrice nx2 dei min-max
  #---------------------------------------------------------
  call.ok<-update.lme.call(my.call, fixed = formulaFix, random=random, data=mf, evaluate=FALSE,
                           control = list(msVerbose = FALSE, niterEM = 100, opt = "optim"))
  if(!is.null(start.pd)) call.ok$random<-quote(list(id=start.pd))
  #--------------------------------------------------------
  kappa0i  <- if(psi.link=="logit") logit(psi.ex,min.Z,max.Z)  else psi.ex #length=n
  if(est.kappa0) kappa0<-mean(kappa0i)
  ki<- kappa0i - kappa0
  etai<- kappa0i + eta.psi
  psi.ex<-if(psi.link=="logit") inv.logit(etai,min.Z,max.Z) else etai  #length=n
  
  #----------------------------------------------------------
  boot.call<-update.lme.call(my.call, y.b~., data=newData, evaluate=FALSE) #salva la call before modifying obj
  it <- 1
  epsilon <- 9
  obj<-o #serve per estrarre la logLik
  b.new<-rep(.1,length(all.vars(formulaFix))) #la risposta conteggiata in all.vars(formulaFix) conta per l'intercetta
  while(abs(epsilon) > tol){
    #if(it==9) browser()
    DD<-if(psi.link=="logit") (max.Z-min.Z)*exp(etai)/((1+exp(etai))^2) else rep(1,n)
    V<-ifelse(Z >psi.ex, -1, 0)
    VD <- V*DD
    mf$U <- pmax(0, Z-psi.ex)
    mf$G0<- rep(delta0i,ni)*VD #rowSums(rep(delta0i,ni)*VD)
    if(id.x.diff){
      Ux<- as.matrix(mf$U*X.diff)
      colnames(Ux)<-nomiUx
      mf[,which(names(mf)%in%nomiUx)]<-Ux
      deltaMatrix<-cbind(rep(delta0i,ni), matrix(delta,nrow=length(V),ncol=length(delta),byrow=TRUE))
      deltaVDx<-deltaMatrix*VD*cbind(1,M.x.diff)
      mf$G0<-rowSums(deltaVDx)
    }
    if(id.z.psi){
      G<-cbind(mf$G0,mf$G0*M.z.psi)
      colnames(G)<-c("G0",nomiG)
      mf[,colnames(G)]<-G
    }
    dev.old <- obj$logLik
    #costruisci l'offset e modifica la risposta..
    Off<- if(est.kappa0)  -kappa0i*mf$G0 else -ki*mf$G0
    if(id.z.psi) Off<- Off - drop(as.matrix(mf[nomiG])%*%kappa[nomiG])
    mf[nomeRispo]<-Rispo-Off
    
    # estimate the model
    ########################################
    obj<-eval(call.ok)
    ########################################
    
    #formulaFix.noG
    #random.noG
    
    b.old<-b.new
    b.new<-fixed.effects(obj)
    ###    if(psi.new>max(Z)| psi.new<min(Z)) stop("estimated psi out of range: try another starting value!")
    dev.new <- obj$logLik#sum((fitted(obj)-my.dd[,paste(formula(obj))[2]])^2) #
    
    
    #===============================================================================
    if (display) {
      flush.console()
      spp <- if (it < 10) " " else NULL
      cat(paste("iter = ", spp, it,
                "  work.LL = ",formatC(dev.new,digits=3,format="f"), #era format="fg"
                "  diff.s = ",formatC(fixef(obj)["U"],digits=3,format="f"), 
                "  kappa0 = ",paste(formatC(fixef(obj)["G0"],digits=3, format="f"), collapse="  "),
                sep=""), "\n")
    }
    
    
    
    #===============================================================================
    epsilon <- abs((dev.new-dev.old)/(dev.old+.1))
    #epsilon <- max(abs((b.new-b.old)/b.old))
    #26/7/16 PERCHE' HO MESSO QUI i CRITERI DI ARRESTO? E' un problema perch? poi il ciclo non
    # termina e i "delta", "kappa0", rimangono quelli dell'iterazione precedente..
    #if(it >= it.max) break
    #if(abs(epsilon) <= tol) break
    it <- it+1
    #stopping rules not met: update the estimates
    ##-------------------------------
    continua<-  (abs(epsilon) > tol && it< it.max)
    #delta0i<-if(inflate.res) inflate.2residuals(obj, coeff=TRUE)[,"U"] else unlist(coef(obj)["U"])    #length=N
    if(id.x.diff) delta <- fixed.effects(obj)[nomiUx]
    
    delta0i <- unlist(coef(obj)["U"])
    kappa0.old <- kappa0 #length=1
    kappa0 <- fixed.effects(obj)["G0"]
    
    if(est.kappa0 && continua){
      kappa0<- if(psi.link=="identity")  adj.psi(kappa0, limZ) else max(min(9,kappa0),-9)
      kappa0 <- kappa0.old + (kappa0 - kappa0.old)*h/2 
      #questo controllo e' sbagliato se link.psi="logit"
      #if(kappa0<= min(Z) || kappa0>=max(Z)) stop("estimated psi outside the range")
    }
    
    
    
    #browser()
    
    kappa0i.old<-kappa0i #length=n
    
    #browser()
    RE<-fn.re(obj) # array n x n.randmEff (2 se sono nested..) x n.termini (U, G0,..) 
    ki<-if("G0" %in% dimnames(RE)[[3]]) rowSums(RE[ , ,"G0", drop=FALSE]) else rep(0,n)
    #NB    RE[ , ,"G0"]  ? una matrice di n.obs righe e che ha in ogni colonna i breakpoint relativi ad ogni livello di nesting.. 
    #      RE[ , J,"G0"] e' l'innermost J=ncol(RE[ , ,"G0"])
    #Quindi i ki sono la somma di tutti i termini random (anche a diversi livelli di nested)
    kappa0i <- kappa0+ki
    
    ########I codici sotto non funzionano con nested r.e.        
    #        ki<-if("G0"%in%names(ranef(obj))) unlist(ranef(obj)["G0"]) else rep(0,N)
    #        kappa0i <- kappa0+ki #length=N
    #        #kappa0i <-if(inflate.res) inflate.2residuals(obj, coeff=TRUE)[,"G0"] else unlist(coef(obj)["G0"]) #length=N
    #        kappa0i<-rep(kappa0i,ni) #+ kappa0i.old #length=n
    #        ki<-rep(ki,ni)
    ###########################
    
    etai<-kappa0i
    if(id.z.psi) {
      kappa.old<-kappa #length=1
      kappa<-fixed.effects(obj)[nomiG]  #esclude G0..
      etai<-etai+drop(Z.psi%*%kappa)
    }
    #questo e' se ci sono parametri con valori *fissati* da non stimare..
    if(anyFixedG){ 
      etai <- etai+ Offs.kappa
    }
    
    #browser()
    
    psi.old <- psi.ex #length=n.obs
    psi.ex<-if(psi.link=="logit") inv.logit(etai,min.Z,max.Z) else etai  #length=n
    #eventuale aggiustamento dei psi.
    #        if(adjust==2){
    #            id.bp<-I(psi.new>minMax[,1]&psi.new<minMax[,2])
    #            psi.new[!id.bp] <- tapply(Z,id,max)[!id.bp]# minMax[!id.bp,2]
    #            }
    
    #if(it==2) browser()
    
    if(it >= (it.max+1)) break
    #        if(abs(epsilon) <= tol) break #NON serve, c'? il while(abs(epsilon) > tol)
    
    #f.pd() la chiamo solo se non ci sono nested r.e. (perch? in quel caso non funziona..) 
    if(J<=1){ #se c'e' SOLO 1 r.e. 
      pd<-f.pd(obj)
      call.ok$random<-quote(list(id=pd))
    }
  } #end_while
  #---------------------------------------------------------------------------------------
  #Adesso devi fare in modo che le linee *veramente si uniscano (no salti), boot restarting e
  #valore di logLik ed infine aggiorna obj<-eval(call.ok)
  
  fixed.noG<-if(is.null(nomiG)) update.formula(call.ok$fixed, paste(".~.-G0",sep="")) 
  else update.formula(call.ok$fixed, paste(".~.-G0-",paste(nomiG, collapse="-"),sep=""))
  if(is.null(random.noG)){ #se "random.noG" non ? stato specificato in segmented.lme()
    random.noG<-random
    #Escludi G0 dalla formula random..
    #  -
    #18/6/16: mi sono reso conto che random pu? essere una lista che contiene diverse formula che includono "G0" (ad es., nel caso di r.e.), quindi "G0" si deve
    # eliminare in ogni formula..
    # Just now I don't know what happen if random is a block matrix.. VERIFICARE.. comunque il codice sotto c'e'..
    
    for(j in 1:J){ #J =n. di random cluster (a des., children %in% school,..)
      #questo sotto ? se random ? una lista e ogni sua componente ha una formula come "attributo".. Dovrebbero rientrare i casi di
      #semplici e nested r.e. NON con una matrice a blocchi..
      if(!is.null(attr(random.noG[[j]], "formula"))){ #semplici e nested r.e.
        if("G0"%in%all.vars(attr(random.noG[[j]], "formula"))){#se la formula della componente j contiene "G0"..
          attr(random.noG[[j]], "formula") <- update.formula(attr(random.noG[[j]], "formula"), ~.-G0)
        }
        #questo sotto e' se c'e' una matrice a blocchi..
      } else {
        #questo sotto e' se c'e' una matrice a blocchi..
        for(k in length(random.noG[[j]])) {
          if(!is.null(attr(random.noG[[j]][[k]], "formula"))){ #Questo ? se ci sono matrici a blocchi quando 
            if("G0"%in%all.vars(attr(random.noG[[j]][[k]], "formula"))){#se la formula della componente j contiene "G0"..
              attr(random.noG[[j]][[k]], "formula") <- update.formula(attr(random.noG[[j]][[k]], "formula"), ~.-G0)
            }
          }
        } #end k=1..K
      }
    } #end j=1..J 
  }
  
  call.ok.noG<-update.lme.call(call.ok, fixed = fixed.noG, random = random.noG)
  mf[nomeRispo]<-Rispo
  obj.noG<-eval(call.ok.noG)
  
  #if(it >= (it.max+1)) warning("max no. of iterations achieved.. refit.boot() suggested", call. = FALSE)
  psi.new<-psi.ex[cumsum(ni)]
  
  #5/7/18: rownames(ranef(obj)[[J]]) sono del tipo "1/1", cio? tengono conto di eventuali nested.. 
  #names(psi.new)<-rownames(ranef(obj)[[J]])
  
  
  #names(psi.new)<-levels(unlist(obj$groups))
  #names(psi.new)<-levels(id)
  ##27/6, nuovo:
  #se id e' numerica levels(id) e' NULL, per cui i psi.new sono senza nomi (e questo da errore in plot.segmented)
  #names(psi.new)<-levels(factor(id)) #funziona anche con nested r.e.??
  #browser()
  rnfGroups<-obj.noG$groups
  
  #names(psi.new)<-levels(rnfGroups[, ncol(rnfGroups)]) #levels ordina per i nomi "nuovi" (se c'? nested 4/10 lo considera prima di 4/9').. 
  names(psi.new)<-rownames(coef(obj.noG)) #oppure unique(rnfGroups[, ncol(rnfGroups)])
  attr(psi.new,which="ni")<-table(rnfGroups[, ncol(rnfGroups)]) 
  
  id.bp<-I(psi.new>=minMax[,1]&psi.new<=minMax[,2])
  attr(psi.new,which="is.break")<-id.bp
  
  #mf$rispo<-Rispo
  #o.new<-lme.formula(rispo ~ x + U + U.x.diff, data = mf, random=list(id=pdDiag(~1+x+U)), method=..)
  #return(o.new)
  
  if(adjust==1){
    #ristima il modello con i nuovi psi ( e le nuove variabili)
    psi.new[!id.bp] <- tapply(Z,id,max)[!id.bp]# minMax[!id.bp,2]
    psi.ex <- rep(psi.new, aa) #length=n.obs
    DD<-fn1(c(rep(kappa0,aa),kappa1), Z.psi ,2, link=psi.link) #length=n.obs
    V<-ifelse(Z >psi.ex, -1, 0)
    my.dd$U<- pmax(0, Z -psi.ex)
    VD <- V*DD
    deltaMatrix<-cbind(rep(betaa,aa), matrix(delta,nrow=length(V),ncol=length(delta),byrow=TRUE))
    deltaVDx<-deltaMatrix*VD*M.x.diff
    G0<-rowSums(deltaVDx)
    G<-G0*M.z.psi
    colnames(G)<-c("G0",paste("G.",colnames(M.z.psi)[-1],collapse="+",sep=""))
    my.dd<-cbind(my.dd, G)
    dev.old <- obj$logLik
    #stima il modello:
    obj<-eval(call.ok)
  }
  
  
  #if(id.z.psi) names(kappa)<- colnames(M.z.psi) #? gi? fatto prima
  RIS <- list("lme.fit"=obj, "lme.fit.noG"=obj.noG, "psi.i"=psi.new, call=match.call())
  if(!is.null(fixed.parms)) RIS$fixed.parms<-fixed.parms
  if(id.z.psi) {
    RIS$fixed.eta.psi<-drop(as.matrix(cbind(1,M.z.psi[cumsum(ni),]))%*%c(kappa0,kappa))
    names(RIS$fixed.eta.psi) <-names(psi.new)
  } else {
    RIS$fixed.eta.psi<-rep(kappa0, length(psi.new))
    names(RIS$fixed.eta.psi) <-names(psi.new)
  }
  if(id.x.diff) {
    RIS$fixed.eta.delta<-drop(as.matrix(cbind(1,M.x.diff[cumsum(ni),]))%*%fixef(obj)[c("U",nomiUx)])
    names(RIS$fixed.eta.delta) <-names(psi.new)
  } else {
    RIS$fixed.eta.delta<- rep(fixef(obj)["U"], length(psi.new))
    names(RIS$fixed.eta.delta) <-names(psi.new)
  }
  
  RIS$fixed.psi<-if(psi.link=="logit") inv.logit(RIS$fixed.eta.psi,min.Z,max.Z) else RIS$fixed.eta.psi
  #browser()
  names(RIS$fixed.psi) <- names(psi.new)
  RIS$call$psi.link<-psi.link #in questo modo il nome e' "completo"..
  RIS$boot.call<-boot.call
  RIS$randomCALL<-randomCALL
  RIS$misc$datacall<- datacall
  #browser()
  #RIS$misc$matrix.psi<- 
  if("G0" %in% dimnames(RE)[[3]]) {
    RIS$misc$matrix.psi<- cbind(fixed=RIS$fixed.psi,drop(RE[cumsum(ni), , "G0", drop = FALSE]))
    colnames(RIS$misc$matrix.psi) <- c("fixed", names(obj$groups))
    rownames(RIS$misc$matrix.psi) <- names(psi.new)#rownames(ranef(obj)[[J]])    
  } else {
    RIS$misc$matrix.psi<- matrix(RIS$fixed.psi, ncol=1) #fixed=RIS$fixed.psi)
    rownames(RIS$misc$matrix.psi) <- names(psi.new)#rownames(ranef(obj)[[J]])    
  }
  
  RIS$namesGZ<-namesGZ
  RIS$Off<-Off
  RIS$rangeZ<- tapply(Z, id, range)
  names(Z)<-id #names(psi.new)
  RIS$Z<-Z
  #browser()
  class(RIS)<- "segmented.lme" #c("segmented.lme","segmented")
  #opz.control<-list(...)
  #if(!is.null(opz.control$n.boot)) n.boot<- opz.control$n.boot
  if(it >= (it.max+1) && n.boot==0) warning("max no. of iterations achieved.. 'n.boot>0' suggested", call. = FALSE)
  if(n.boot>0){
    if(display) cat("Implementing bootstrap restarting..\n")
    RIS <- reboot.slme(RIS, B=n.boot, display=display, break.boot=control$break.boot ,
                       seed=control$seed, msg=display)#, metodo=1, frac=1, it.max=6, it.max.b=5, start=NULL, msg=TRUE)
  }
  RIS
}

