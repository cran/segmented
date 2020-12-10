`segmented.lm` <-
function(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(), model = TRUE, keep.class=FALSE, ...) {
  
  build.all.psi<-function(psi, fixed.psi){
    all.names.psi<-union(names(psi),names(fixed.psi))
    all.psi<-vector("list", length=length(all.names.psi))
    names(all.psi)<- all.names.psi
    for(i in names(all.psi)) {
      if(!is.null(psi[[i]])){
        psi[[i]]<-sort(psi[[i]])
        names(psi[[i]])<-paste("U",1:length(psi[[i]]),".",i,sep="")
      }
      if(!is.null(fixed.psi[[i]])){
        fixed.psi[[i]]<-sort(fixed.psi[[i]])
        names(fixed.psi[[i]])<-	paste("U",1:length(fixed.psi[[i]]),".fixed.",i,sep="")
      }
      all.psi[[i]]<-sort(c(psi[[i]],fixed.psi[[i]]))
    }
    return(all.psi)
  }
  ##===inizio funzione
  
  if(missing(seg.Z)) {
      if(length(all.vars(formula(obj)))==2) seg.Z<- as.formula(paste("~", all.vars(formula(obj))[2])) else stop("please specify 'seg.Z'")
  }
  n.Seg<-length(all.vars(seg.Z))
  id.npsi<-FALSE
  
  if("V" %in% sub("V[1-9]*[0-9]","V", c(all.vars(seg.Z), all.vars(formula(obj) )[-1]))) stop("variable names 'V', 'V1', .. are not allowed")
  if("U" %in% sub("U[1-9]*[0-9]","U", c(all.vars(seg.Z), all.vars(formula(obj) )[-1]))) stop("variable names 'U', 'U1', .. are not allowed")
  if(any(c("$","[") %in% all.names(seg.Z))) stop(" '$' or '[' not allowed in 'seg.Z' ")
  
    if(missing(psi)) { 
     if(n.Seg==1){
       if(missing(npsi)) npsi<-1
       npsi<-lapply(npsi, function(.x).x)
       if(length(npsi)!=length(all.vars(seg.Z))) stop("seg.Z and npsi do not match") 
       names(npsi)<-all.vars(seg.Z)
     } else {#se n.Seg>1
       #if(missing(npsi)) stop(" with multiple segmented variables in seg.Z, 'psi' or 'npsi' should be supplied", call.=FALSE) 
       if (missing(npsi)) {
         npsi<-rep(1, n.Seg)
         names(npsi)<-all.vars(seg.Z)
       }
       if(length(npsi)!=n.Seg) stop(" 'npsi' and seg.Z should have the same length")
       if(!all(names(npsi) %in% all.vars(seg.Z))) stop(" names in 'npsi' and 'seg.Z' do not match")    
       }
     psi<-lapply(npsi, function(.x) rep(NA,.x))
     id.npsi<-TRUE ##id.npsi<-FALSE #e' stato fornito npsi?
    } else {
     if(n.Seg==1){
        if(!is.list(psi)) {psi<-list(psi);names(psi)<-all.vars(seg.Z)}
        } else {#se n.Seg>1
        if(!is.list(psi)) stop("with multiple terms in `seg.Z', `psi' should be a named list")
        if(n.Seg!=length(psi)) stop("A wrong number of terms in `seg.Z' or `psi'")
        if(!all(names(psi)%in%all.vars(seg.Z))) stop("Names in `seg.Z' and `psi' do not match")
        }
    }


    fc<- min(max(abs(control$fc),.8),1)       
    min.step<-control$min.step
    alpha<-control$alpha
    it.max <- old.it.max<- control$it.max
    digits<-control$digits
    toll <- control$toll
    if(toll<0) stop("Negative tolerance ('tol' in seg.control()) is meaningless", call. = FALSE)
    visual <- control$visual
    
    stop.if.error<-control$stop.if.error
    fix.npsi<-fix.npsi<-control$fix.npsi
    if(!is.null(stop.if.error)) {#if the old "stop.if.error" has been used..
      warning(" Argument 'stop.if.error' is working, but will be removed in the next releases. Please use 'fix.npsi' for the future..")
    } else {
      stop.if.error<-fix.npsi
    }
    
    n.boot<-control$n.boot
    size.boot<-control$size.boot
    gap<-control$gap
    random<-control$random
    pow<-control$pow
    conv.psi<-control$conv.psi
    visualBoot<-FALSE
    if(n.boot>0){
        if(!is.null(control$seed)) {
            set.seed(control$seed)
            employed.Random.seed<-control$seed
              } else {
            employed.Random.seed<-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
            set.seed(employed.Random.seed)
              }
        if(visual) {visual<-FALSE; visualBoot<-TRUE}# warning("`display' set to FALSE with bootstrap restart", call.=FALSE)}
#        if(!stop.if.error) stop("Bootstrap restart only with a fixed number of breakpoints")
     }
    last <- control$last
    K<-control$K
    h<-control$h
    #-------------------------------
    orig.call<-Call<-mf<-obj$call
    orig.call$formula<- mf$formula<-formula(obj) #per consentire lm(y~.)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if(class(mf$formula)=="name" && !"~"%in%paste(mf$formula)) mf$formula<-eval(mf$formula)
    #orig.call$formula<-update.formula(orig.call$formula, paste("~.-",all.vars(seg.Z)))
#
#genn 2013. dalla versione 0.2.9-4 ho tolto if(length(.. Tra l'altro non capisco perche' lo avevo fatto
#    if(length(all.vars(formula(obj)))>1){
#    mf$formula<-update.formula(mf$formula,paste(paste(seg.Z,collapse=".+"),"+",paste(all.vars(formula(obj))[-1],collapse="+")))
#    } else {
#    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
#    }
#nov 2013 dalla versione 0.3-0.0 (che dovrebbe essere successiva alla 0.2-9.5) viene creato anche il modelframe esteso che comprende
# termini "originali", prima che fossero trasformati (Ad es., x prima che ns(x) costruisca le basi). Questo permette di avere termini
# ns(), poly(), bs() nel modello di partenza
    mfExt<- mf
    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
    #mfExt$formula<- update.formula(mfExt$formula,paste(paste(seg.Z,collapse=".+"),"+",paste(all.vars(formula(obj)),collapse="+")))   
#    mfExt$formula<- if(!is.null(obj$call$data)) 
#        update.formula(mf$formula,paste(".~",paste(all.vars(obj$call), collapse="+"),"-",obj$call$data,sep=""))
#            else update.formula(mf$formula,paste(".~",paste(all.vars(obj$call), collapse="+"),sep=""))
#-----------
#   browser()
   if(!is.null(obj$call$offset) || !is.null(obj$call$weights) || !is.null(obj$call$subset)){ 
      mfExt$formula <- 
          update.formula(mf$formula, 
          paste(".~.+", paste(
          c(all.vars(obj$call$offset), 
            all.vars(obj$call$weights),
            all.vars(obj$call$subset)), 
            collapse = "+")
            ))
          }
    mf <-  eval(mf, parent.frame())
    n<-nrow(mf)
    #questo serve per inserire in mfExt le eventuali variabili contenute nella formula con offset(..)
    nomiOff<-setdiff(all.vars(formula(obj)), names(mf))
    if(length(nomiOff)>=1) mfExt$formula<-update.formula(mfExt$formula,paste(".~.+", paste( nomiOff, collapse="+"), sep=""))
#----------------------------------------------------
#ago 2014 c'e' la questione di variabili aggiuntive...
nomiTUTTI<-all.vars(mfExt$formula) #comprende anche altri nomi (ad es., threshold) "variabili"
nomiNO<-NULL 
for(i in nomiTUTTI){
    r<-try(eval(parse(text=i), parent.frame()), silent=TRUE)
    if(class(r)!="try-error" && length(r)==1 && !is.function(r) && !i%in%names(mf)) nomiNO[[length(nomiNO)+1]]<-i
    }
#nomiNO dovrebbe contenere i nomi delle "altre variabili" (come th in subset=x<th) 
if(!is.null(nomiNO)) mfExt$formula<-update.formula(mfExt$formula,paste(".~.-", paste( nomiNO, collapse="-"), sep=""))
#funziona ma se c'e' qualche variabile con il nome di una funzione (ad es., "filter")
# anche questa variabile va a finire in nomiNO e quindi poi viene eliminata
# problema risolto con "r<-try(evalq(.."? (prima era "r<-try(eval(..")
#----------------------------------------------------
    
    #nomiAgg<- setdiff(names(mf), names(mfExt))
    #if(length(nomiAgg)>0) mfExt$formula<-update.formula(mfExt$formula,paste(".~.+", paste( nomiAgg, collapse="+"), sep=""))
    
    mfExt<-eval(mfExt, parent.frame())
    
    #mantieni in mfExt solo le variabili che NON ci sono in mf (cosi la funzione occupa meno spazio..)
    #mfExt<-mfExt[,setdiff(names(mfExt), names(mf)),drop=FALSE]
    
    weights <- as.vector(model.weights(mf))
    offs <- as.vector(model.offset(mf))
    mt <- attr(mf, "terms")
    interc<-attr(mt,"intercept")
    y <- model.response(mf, "any")
    XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf, obj$contrasts)

    #il cambio in mf da "offset(_nomevar_)" al "_nomevar_" deve avvenire dopo "model.matrix(mt, mf, contrasts)" 
#    if(!is.null(offs)){
#      #id.offs<-pmatch("offset",names(mf)) #questa identifica il nome offset(..). ELiminarlo dal dataframe? non conviene altrimenti nel model.frame non risulta l'offset
#      id.offs<- which(grepl("(offset)", names(mf))) #per consentire anche offset come argomento di glm()
#      names(mf)[id.offs]<- all.vars(formula(paste("~", names(mf)[id.offs])), functions=FALSE)
#      }

    namesXREG0<-colnames(XREG)
    #nameLeftSlopeZero<-setdiff(all.vars(seg.Z), all.vars(formula(obj)))
    nameLeftSlopeZero<-setdiff(all.vars(seg.Z), names(coef(obj))) #in questo modo riconosce che sin(x*pi) NON e' x, ad esempio.
    namesXREG0<-setdiff(namesXREG0, nameLeftSlopeZero)
    id.duplic<-match(all.vars(formula(obj)),all.vars(seg.Z),nomatch=0)>0
    if(any(id.duplic)) {
        #new.mf<-mf[,id.duplic,drop=FALSE]
        new.mf<-mf[,all.vars(formula(obj))[id.duplic],drop=FALSE]
        new.XREGseg<-data.matrix(new.mf)
        XREG<-cbind(XREG,new.XREGseg)
        }
    n.psi<- length(unlist(psi))
    id.n.Seg<-(ncol(XREG)-n.Seg+1):ncol(XREG)
    XREGseg<-XREG[,id.n.Seg,drop=FALSE]
    #XREG<-XREG[,-id.n.Seg,drop=FALSE]
    #XREG<-model.matrix(obj0) non va bene perche' non elimina gli eventuali mancanti in seg.Z..
    #Due soluzioni
    #XREG<-XREG[,colnames(model.matrix(obj)),drop=FALSE]
    #XREG<-XREG[,match(c("(Intercept)",all.vars(formula(obj))[-1]),colnames(XREG),nomatch =0),drop=FALSE]
    XREG <- XREG[, match(c("(Intercept)", namesXREG0),colnames(XREG), nomatch = 0), drop = FALSE]
    XREG<-XREG[,unique(colnames(XREG)), drop=FALSE]
    #################
    #if(ncol(XREGseg)==1 && length(psi)==1 && n.psi==1 && !any(is.na(psi))) { if(psi==Inf) psi<-median(XREGseg)}
    #################
    n <- nrow(XREG)
    Z<-lapply(apply(XREGseg,2,list),unlist) #prende anche i nomi!
    name.Z <- names(Z) <- colnames(XREGseg)
    if(length(Z)==1 && is.vector(psi) && (is.numeric(psi)||is.na(psi))){
        psi <- list(as.numeric(psi))
        names(psi)<-name.Z
        }
    if (!is.list(Z) || !is.list(psi) || is.null(names(Z)) || is.null(names(psi))) stop("Z and psi have to be *named* list")
    id.nomiZpsi <- match(names(Z), names(psi))
    if ((length(Z)!=length(psi)) || any(is.na(id.nomiZpsi))) stop("Length or names of Z and psi do not match")
    nome <- names(psi)[id.nomiZpsi]
    psi <- psi[nome]
    if(id.npsi){
      for(i in 1:length(psi)) {
        K<-length(psi[[i]])
        if(any(is.na(psi[[i]]))) psi[[i]]<-if(control$quant) {quantile(Z[[i]], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)} else {(min(Z[[i]])+ diff(range(Z[[i]]))*(1:K)/(K+1))}
      }
    } else {
    for(i in 1:length(psi)) {
        if(any(is.na(psi[[i]]))) psi[[i]]<-if(control$quant) {quantile(Z[[i]], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)} else {(min(Z[[i]])+ diff(range(Z[[i]]))*(1:K)/(K+1))}
        }
    }

  #########==================== SE PSI FIXED
  id.psi.fixed <- FALSE
  if(!is.null(fixed.psi)){
    id.psi.fixed <- TRUE
    if(is.numeric(fixed.psi) && n.Seg==1) {
      fixed.psi<-list(fixed.psi)
      names(fixed.psi)<-all.vars(seg.Z)
    }
    if(is.list(fixed.psi)) {
      if(!(names(fixed.psi) %in% all.vars(seg.Z))) stop("names(fixed.psi) is not a subset of variables in 'seg.Z' ")
    } else {
      stop(" 'fixed.psi' has to be a named list ")
      } 
    fixed.psi<-lapply(fixed.psi, sort)
    Zfixed<-matrix(unlist(mapply(function(x,y)rep(x,y),Z[names(fixed.psi)], sapply(fixed.psi, length), SIMPLIFY = TRUE)), nrow=n)
    n.fixed.psi<-sapply(fixed.psi, length)
    rip.nomi <- rep( names(fixed.psi), n.fixed.psi)
    rip.numeri <- unlist(lapply(n.fixed.psi, function(.x) 1:.x))
    colnames(Zfixed) <- paste("U", rip.numeri,".fixed.",rip.nomi, sep="")
    PSI <- matrix(unlist(fixed.psi), ncol=ncol(Zfixed), nrow=n, byrow = TRUE)
    fixedU<-(Zfixed-PSI)*(Zfixed>PSI)
    XREG<-cbind(XREG, fixedU)
  }
#########====================END  SE PSI FIXED

    initial.psi<-psi
    a <- sapply(psi, length)
    #per evitare che durante il processo iterativo i psi non siano ordinati
    id.psi.group <- rep(1:length(a), times = a) #identificativo di apparteneza alla variabile
    #

    #Znew <- list()
    #for (i in 1:length(psi)) Znew[[length(Znew) + 1]] <- rep(Z[i], a[i])
    #Z <- matrix(unlist(Znew), nrow = n)
    Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n)
    psi <- unlist(psi)
    #se psi e' numerico, la seguente linea restituisce i valori ordinati all'interno della variabile..
    psi<-unlist(tapply(psi,id.psi.group,sort))
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
    #controllo se psi e' ammissibile..
    c1 <- apply((Z <= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo <)
    c2 <- apply((Z >= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo >)
    if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2)) ) stop("starting psi out of the admissible range")
    
    colnames(Z) <- nomiZ <- rep(nome, times = a)
    ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))], function(.x) {1:.x})))
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ, sep = ".")
    nomiV <- paste("V", ripetizioni, sep = "")
    nomiV <- paste(nomiV, nomiZ, sep = ".")
    #forse non serve crearsi l'ambiente KK, usa mf..
    #obj <- update(obj, formula = Fo, data = mf)
    #if (model.frame) obj$model <- mf
  #controlla che model.frame() funzioni sull'oggetto restituito
#    KK <- new.env()
#    for (i in 1:ncol(objframe$model)) assign(names(objframe$model[i]), objframe$model[[i]], envir = KK)
    if (it.max == 0) {
        #mf<-cbind(mf, mfExt)
        U <- (Z>PSI)*(Z-PSI) #pmax((Z - PSI), 0)
        colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
        nomiU <- paste("U", colnames(U), sep = "")
        #for (i in 1:ncol(U)) assign(nomiU[i], U[, i], envir = KK)
        #e' necessario il for? puoi usare colnames(U)<-nomiU;mf[nomiU]<-U
        for(i in 1:ncol(U)) mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
        Fo <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
        obj <- update(obj, formula = Fo, evaluate=FALSE, data=mfExt) #data = mf, 
        if(!is.null(obj[["subset"]])) obj[["subset"]]<-NULL
        obj<-eval(obj, envir=mfExt)
        if (model) obj$model <-mf  #obj$model <- data.frame(as.list(KK))

        psi <- cbind(psi, psi, 0)
        rownames(psi) <- paste(paste("psi", ripetizioni, sep = ""), nomiZ, sep=".")
        colnames(psi) <- c("Initial", "Est.", "St.Err")

        #names(psi)<-paste(paste("psi", ripetizioni, sep = ""), nomiZ, sep=".")
        obj$psi <- psi
        return(obj)
    }
    #XREG <- model.matrix(obj) creata sopra
    #o <- model.offset(objframe)
    #w <- model.weights(objframe)
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(offs)) offs <- rep(0, n)
    initial <- psi
    obj0 <- obj
    dev0<-sum(obj$residuals^2)
    list.obj <- list(obj)
#    psi.values <- NULL
    nomiOK<-nomiU
    invXtX<-chol2inv(qr.R(obj$qr)) #(XtX)^{-1}
    Xty<-crossprod(XREG,y)
    opz<-list(toll=toll,h=h, stop.if.error=stop.if.error, dev0=dev0, visual=visual, it.max=it.max,
        nomiOK=nomiOK, id.psi.group=id.psi.group, gap=gap, visualBoot=visualBoot, pow=pow, digits=digits,invXtX=invXtX, Xty=Xty, 
        conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step, fc=fc)
    if(n.boot<=0){
    obj<-seg.lm.fit(y,XREG,Z,PSI,weights,offs,opz)
    } else {
    obj<-seg.lm.fit.boot(y, XREG, Z, PSI, weights, offs, opz, n.boot=n.boot, size.boot=size.boot, random=random) #jt, nonParam
      }
    if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(obj0)
        }
    
    if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
    id.psi.group<-obj$id.psi.group
    nomiOK<-obj$nomiOK
    #nomiFINALI<-unique(sapply(strsplit(nomiOK, split="[.]"), function(x)x[2])) #nomi delle variabili con breakpoint stimati!
    #nomiFINALI<-sub("U[1-9].", "", nomiOK) #nomi originali delle variabili con breakpoint stimati!
    nomiFINALI<- unique(sub("U[1-9]*[0-9].", "", nomiOK))
    #se e' stata usata una proc automatica "nomiFINALI" sara' differente da "name.Z"
    nomiSenzaPSI<-setdiff(name.Z,nomiFINALI)
    if(length(nomiSenzaPSI)>=1) warning("no breakpoints found for: ", paste(nomiSenzaPSI," "), call. = FALSE)
    it<-obj$it
    psi<-obj$psi
    psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
    U<-obj$U
    V<-obj$V
    id.warn<-obj$id.warn
    rangeZ<-obj$rangeZ
    obj<-obj$obj
    k<-length(psi)
    beta.c<-coef(obj)[paste("U", 1:ncol(U), sep = "")]
    #psi.values[[length(psi.values) + 1]] <- psi #non c'e' bisogno!
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))

    #se usi una procedura automatica devi cambiare ripetizioni, nomiU e nomiV, e quindi:
    length.psi<-tapply(as.numeric(as.character(names(psi))), as.numeric(as.character(names(psi))), length)
    forma.nomiU<-function(xx,yy)paste("U",1:xx, ".", yy, sep="")
    forma.nomiVxb<-function(xx,yy)paste("psi",1:xx, ".", yy, sep="")
    
    #nomiU   <- unlist(mapply(forma.nomiU, length.psi, name.Z)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
    #nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, name.Z))
    nomiU   <- unlist(mapply(forma.nomiU, length.psi, nomiFINALI)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
    nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, nomiFINALI))

    #########========================= SE PSI FIXED
    psi.list<-vector("list", length=length(unique(nomiZ)))
    names(psi.list)<-unique(nomiZ)
    #names(psi)<-nomiZ #se e' una procedure automatica nomiZ puo essere piu lungo dei breakpoints "rimasti" 
    names(psi)<-rep(nomiFINALI, length.psi)
    for(i in names(psi.list)){
      psi.list[[i]]<-psi[names(psi)==i]
    }

    #mf<-cbind(mf, mfExt) #questo creava ripetizioni..
    for(i in 1:ncol(U)) {
        mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
        mfExt[nomiVxb[i]]<-mf[nomiVxb[i]]<-Vxb[,i]
        }
    nnomi <- c(nomiU, nomiVxb)
#    browser()
    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
    #########========================= SE PSI FIXED
    if(id.psi.fixed){
      for(i in 1:ncol(fixedU)) mfExt[colnames(fixedU)[i]]<-mf[colnames(fixedU)[i]]<-fixedU[,i]
      Fo<-update.formula(Fo, paste(c("~.",colnames(fixedU)), collapse="+"))
    }
    #objF <- update(obj0, formula = Fo, data = KK)
    objF <- update(obj0, formula = Fo,  evaluate=FALSE, data = mfExt)
    #eliminiamo subset, perche' se e' del tipo subset=x>min(x) allora continuerebbe a togliere 1 osservazione 
    if(!is.null(objF[["subset"]])) objF[["subset"]]<-NULL
    objF<-eval(objF, envir=mfExt)
    objF$offset<- obj0$offset
    isNAcoef<-any(is.na(objF$coefficients))
    if(isNAcoef){
      if(stop.if.error) {
       cat("breakpoint estimate(s):", as.vector(psi),"\n")
       stop("at least one coef is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)", 
         call. = FALSE)
          } else {
        warning("some estimate is NA: too many breakpoints? 'var(hat.psi)' cannot be computed \n ..returning a 'lm' model", call. = FALSE)
        Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
        objF <- update(obj0, formula = Fo,  evaluate=TRUE, data = mfExt)
        names(psi)<-nomiVxb
        objF$psi<-psi
        return(objF)      
        }
    }
    if(!gap){
      names.coef<-names(objF$coefficients)
      #questi codici funzionano e si basano sull'assunzioni che le U e le V siano ordinate..
      names(obj$coefficients)[match(c(paste("U",1:k, sep=""), paste("V",1:k, sep="")), names(coef(obj)))]<- nnomi  
      objF$coefficients[names.coef]<-obj$coefficients[names.coef] #sostituisce gli 0 
      #objF$coefficients<-obj$coefficients
      #names(objF$coefficients)<-names.coef
      objF$fitted.values<-obj$fitted.values
      objF$residuals<-obj$residuals
      #objF$qr<-obj$qr #NON credo..
    }
    Cov <- vcov(objF) 
    id <- match(nomiVxb, names(coef(objF)))
    vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
    #if(length(initial)!=length(psi)) initial<-rep(NA,length(psi))
    a<-tapply(id.psi.group, id.psi.group, length) #ho sovrascritto "a" di sopra, ma non dovrebbe servire..
    
    ris.psi<-matrix(NA,length(psi),3)
    colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
    rownames(ris.psi) <- nomiVxb
    ris.psi[,2]<-psi
    ris.psi[,3]<-sqrt(vv)
#NB "a" deve essere un vettore che si appatta con "initial.psi" per ottnetere "initial" sotto... Se una variabile alla fine risulta
# senza breakpoint questo non avviene e ci sono problemi nella formazione di "initial". Allora costruisco a.ok
    a.ok<-NULL
    for(j in name.Z){
        if(j %in% nomiFINALI) {
          a.ok[length(a.ok)+1]<-a[1]
          a<-a[-1]
          } else {
          a.ok[length(a.ok)+1]<-0
          } #ifelse(name.Z %in% nomiFINALI,1,0)
        }
    #initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi, a.ok, SIMPLIFY = TRUE))
    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi[nomiFINALI], a.ok[a.ok!=0], SIMPLIFY = TRUE))
    if(stop.if.error)  ris.psi[,1]<-initial 
    #psi <- cbind(initial, psi, sqrt(vv))
    #rownames(psi) <- colnames(Cov)[id]
    
    objF$rangeZ <- rangeZ
    objF$psi.history <- psi.values
    objF$psi <- ris.psi
    objF$it <- it 
    objF$epsilon <- obj$epsilon
    objF$call <- match.call()
    objF$nameUV <- list(U = drop(nomiU), V = rownames(ris.psi), Z = nomiFINALI) #Z = name.Z
    objF$id.group <- if(length(name.Z)<=1) -rowSums(as.matrix(V))
    objF$id.psi.group <- id.psi.group
    objF$id.warn <- id.warn
    objF$orig.call<-orig.call
    objF$indexU<-build.all.psi(psi.list, fixed.psi)
    if(model)  objF$model <- mf #objF$mframe <- data.frame(as.list(KK))
    if(n.boot>0) objF$seed<-employed.Random.seed
    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last)
        list.obj <- list.obj[[length(list.obj)]]
    return(list.obj)
}
