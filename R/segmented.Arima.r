segmented.Arima<-
function(obj, seg.Z, psi, npsi, control = seg.control(), model = TRUE, keep.class=FALSE, ...) {
#Richiede control$f.obj that should be a string like "sum(x$residuals^2)" or "x$dev"
#-----------------
  dpmax<-function(x,y,pow=1){
    #deriv pmax
    if(pow==1) -(x>y) #ifelse(x>y, -1, 0)
    else -pow*((x-y)*(x>y))^(pow-1)#-pow*pmax(x-y,0)^(pow-1)
  }
#-----------
    # n.Seg<-1
    # if(missing(seg.Z) && length(all.vars(o$call$xreg))==1) seg.Z<- as.formula(paste("~", all.vars(o$call$xreg)))
    # if(missing(psi)){if(length(all.vars(seg.Z))>1) stop("provide psi") else psi<-Inf}
    # if(length(all.vars(seg.Z))>1 & !is.list(psi)) stop("`psi' should be a list with more than one covariate in `seg.Z'")
    # if(is.list(psi)){
    #   if(length(all.vars(seg.Z))!=length(psi)) stop("A wrong number of terms in `seg.Z' or `psi'")
    #   if(any(is.na(match(all.vars(seg.Z),names(psi), nomatch = NA)))) stop("Variables in `seg.Z' and `psi' do not match")
    #   n.Seg <- length(psi)
    #   }
    # if(length(all.vars(seg.Z))!=n.Seg) stop("A wrong number of terms in `seg.Z' or `psi'")
  if(missing(seg.Z)) {
    if(length(all.vars(formula(obj)))==2) seg.Z<- as.formula(paste("~", all.vars(formula(obj))[2])) else stop("please specify 'seg.Z'")
  }
  if("V" %in% sub("V[1-9]*[0-9]","V", c(all.vars(seg.Z), names(coef(obj))))) stop("variable names 'V', 'V1', .. are not allowed")
  if("U" %in% sub("U[1-9]*[0-9]","U", c(all.vars(seg.Z), names(coef(obj))))) stop("variable names 'U', 'U1', .. are not allowed")
  if(any(c("$","[") %in% all.names(seg.Z))) stop(" '$' or '[' not allowed in 'seg.Z' ")
  
  n.Seg<-length(all.vars(seg.Z))
  id.npsi<-FALSE
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
      if(!stop.if.error) stop("Bootstrap restart only with a fixed number of breakpoints")
    }
    last <- control$last
    K<-control$K
    h<-control$h
    #    if(h<1) it.max<-it.max+round(it.max/2)
    name.Z <-all.vars(seg.Z)
    if(length(name.Z)!=n.Seg) stop("errore strano 1")
    Z<-sapply(name.Z, function(xx) eval(parse(text=xx))) #e' sempre una matrice
    if(length(name.Z)!=ncol(Z)) stop("errore strano 2")
    n<-nrow(Z)
    n.psi<- length(unlist(psi))
    
    #################
    #if(ncol(Z)==1 && length(psi)==1 && n.psi==1 && !any(is.na(psi))) { if(psi==Inf) psi<-median(Z)}
    #################
    
    if(ncol(Z)==1 && is.vector(psi) && (is.numeric(psi)||is.na(psi))){
      psi <- list(as.numeric(psi))
      names(psi)<-name.Z
    }
    if (!is.list(psi) || is.null(names(psi))) stop("psi should be a *named* list")
    id.nomiZpsi <- match(colnames(Z), names(psi))
    if ((ncol(Z)!=length(psi)) || any(is.na(id.nomiZpsi))) stop("Length or names of Z and psi do not match")
    nome <- names(psi)[id.nomiZpsi]
    psi <- psi[nome]
    if(id.npsi){
      for(i in 1:length(psi)) {
        K<-length(psi[[i]])
        if(any(is.na(psi[[i]]))) psi[[i]]<-if(control$quant) {quantile(Z[,i], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)} else {(min(Z[,i])+ diff(range(Z[,i]))*(1:K)/(K+1))}
      }
    } else {
      for(i in 1:length(psi)) {
        if(any(is.na(psi[[i]]))) psi[[i]]<-if(control$quant) {quantile(Z[,i], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)} else {(min(Z[,i])+ diff(range(Z[,i]))*(1:K)/(K+1))}
      }
    }
    initial.psi<-psi
    a <- sapply(psi, length)
    #per evitare che durante il processo iterativo i psi non siano ordinati
    id.psi.group <- rep(1:length(a), times = a) #identificativo di apparteneza alla variabile
    
    Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n,byrow = TRUE)
    #negli altri metodi Z e' una lista per cui la linea di sopra diventa
    #Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n)
    colnames(Z) <- nomiZ.vett <- rep(nome, times = a) #SERVE??? si perche' Z e' senza colnames
    
    psi <- unlist(psi)
    #se psi e' numerico, la seguente linea restituisce i valori ordinati all'interno della variabile..
    psi<-unlist(tapply(psi,id.psi.group,sort))
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)

        
    #controllo se psi e' ammissibile..
    c1 <- apply((Z <= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo <)
    c2 <- apply((Z >= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo >)
    if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2)) ) stop("starting psi out of the admissible range")
    
    #ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))], function(xxx) {1:xxx})))
    ripetizioni <- as.vector(unlist(tapply(id.psi.group, id.psi.group, function(x) 1:length(x) )))
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ.vett, sep = ".")
    nomiV <- paste("V", ripetizioni, sep = "")
    nomiV <- paste(nomiV, nomiZ.vett, sep = ".")
    nnomi <- c(nomiU, nomiV)
    
    XREG<-eval(obj$call$xreg)
    if(!is.null(XREG)){
      #se ci sono factor?
      nomiXREG<-setdiff(names(obj$coef),c("intercept", paste("ar",1:100,sep=""), paste("ma",1:100,sep=""), 
                                          paste("sma",1:100,sep=""), paste("sar",1:100,sep="")))
      XREG<-matrix(XREG, ncol=length(nomiXREG))
      colnames(XREG)<-nomiXREG
      #if((""%in%colnames(XREG)) || (" "%in%colnames(XREG))) stop("all columns in the matrix 'xreg' of 'obj' should be named.. ")
      if(length(nomiXREG) != ncol(XREG)) stop("ncol(XREG) does not match names of regression coefficients")
    }
    mio.init<-mio.init.noV<-NULL
    X<-NULL
    call.ok <- update(obj,  xreg = X, init=mio.init, evaluate=FALSE) #ho messo X, piuttosto che cbind(XREG,U,V) 
    call.noV <- update(obj, xreg = cbind(XREG,U), init=mio.init.noV,  evaluate=FALSE) #, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
    #    call.noV <- update(obj, formula = Fo.noV,  evaluate=FALSE, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
    
    if (it.max == 0) {
      U<-(Z-PSI)*(Z>PSI)
      colnames(U)<-nomiU
      obj1 <- eval(call.noV) #, envir=mfExt)
      return(obj1)
    }
    
    #obj1 <- eval(call.ok, envir=mfExt)
    initial <- psi
    obj0 <- obj
    
    dev0<- -obj$loglik
    if(is.na(dev0)) dev0<-10
    
    list.obj <- list(obj)
    nomiOK<-nomiU
    
    opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,dev0=dev0,visual=visual,it.max=it.max,
              nomiOK=nomiOK, id.psi.group=id.psi.group, gap=gap, visualBoot=visualBoot, pow=pow, digits=digits,
              conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step,fc=fc)
    
    opz$call.ok<-call.ok
    opz$call.noV<-call.noV
    opz$nomiU<-nomiU
    opz$nomiV<-nomiV

    if(n.boot<=0){
      obj<-seg.Ar.fit(obj, XREG, Z, PSI, opz)
    } else {
      obj<-seg.Ar.fit.boot(obj, XREG, Z, PSI, opz, n.boot=n.boot, size.boot=size.boot, random=random) #jt, nonParam
    }
    
    if(!is.list(obj)){
      warning("No breakpoint estimated", call. = FALSE)
      return(obj0)
    }
    id.warn<-obj$id.warn
    id.psi.group<-obj$id.psi.group
    nomiU<-nomiOK<-obj$nomiOK #sarebbe nomiU
    #--
    nomiVxb<-sub("U","psi", nomiOK) #nomiVxb<-paste("psi",sapply(strsplit(nomiOK,"U"), function(x){x[2]}), sep="")
    nomiFINALI<-unique(sub("U[1-9]*[0-9].", "", nomiOK))     #nomiFINALI<-unique(sapply(strsplit(nomiOK, split="[.]"), function(x)x[2])) #nomi delle variabili con breakpoint stimati!
    #se e' stata usata una proc automatica "nomiFINALI" sara' differente da "name.Z"
    nomiSenzaPSI<-setdiff(name.Z,nomiFINALI)
    if(length(nomiSenzaPSI)>=1) warning("no breakpoints found for: ", paste(nomiSenzaPSI," "), call. = FALSE)
    
    #--
    it<-obj$it
    psi<-obj$psi
    psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
    U<-obj$U
    V<-obj$V
    #    return(obj)
    
    #if(any(table(rowSums(V))<=1)) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close")
    for(jj in colnames(V)) {
      VV<-V[, which(colnames(V)==jj), drop=FALSE]
      sumV<-abs(rowSums(VV))
      if( #(any(diff(sumV)>=2)|| #se ci sono due breakpoints uguali
        any(table(sumV)<=1) && stop.if.error) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
    }
    rangeZ<-obj$rangeZ
    obj<-obj$obj
    k<-length(psi)
    all.coef<-coef(obj)
    names(all.coef)<-c(names(obj0$coef), nomiU, nomiVxb)
    beta.c<- all.coef[nomiU]
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
    nnomi <- c(nomiU, nomiVxb)
    XREG.ok<-cbind(XREG, U, Vxb)
    colnames(XREG.ok)[((ncol(XREG.ok)-length(nnomi)+1):ncol(XREG.ok))]<- nnomi
    objF <- update(obj0,  xreg = XREG.ok, evaluate=TRUE) 

    #    #se usi una procedura automatica devi cambiare ripetizioni, nomiU e nomiV, e quindi:
    #    length.psi<-tapply(as.numeric(as.character(names(psi))), as.numeric(as.character(names(psi))), length)

    if(any(is.na(objF$coef)) && stop.if.error){
      stop("at least one coef estimate is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)", call. = FALSE)
    }
      names.coef <- names(coef(objF))
      #names(obj$coef)<- names.coef# all.coef ha gia' i nomi..
      objF$coef[names.coef]<-all.coef[names.coef]
      objF$residuals<- obj$residuals
      objF$loglik<-obj$loglik
      objF$sigma2 <-obj$sigma2
      objF$aic <- obj$aic + 2*k
      
      if(any(is.na(objF$coef))){ 
        stop("some estimate is NA: premature stopping with a large number of breakpoints?", call. = FALSE)
      }
    
    Cov<-objF$var.coef
    vv<- Cov[nomiVxb, nomiVxb, drop=FALSE]
    ris.psi<-matrix(NA,length(psi),3)
    colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
    rownames(ris.psi) <- nomiVxb
    ris.psi[,2]<-psi
    ris.psi[,3]<-sqrt(diag(vv))
    a<-tapply(id.psi.group, id.psi.group, length) #ho sovrascritto "a" di sopra, ma non dovrebbe servire..
    a.ok<-NULL
    for(j in name.Z){
      if(j %in% nomiFINALI) {
        a.ok[length(a.ok)+1]<-a[1]
        a<-a[-1]
      } else {
        a.ok[length(a.ok)+1]<-0
      } #ifelse(name.Z %in% nomiFINALI,1,0)
    }
    #    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi, a.ok, SIMPLIFY = TRUE))
    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi[nomiFINALI], a.ok[a.ok!=0], SIMPLIFY = TRUE))
    if(stop.if.error)  ris.psi[,1]<-initial
    
    objF$Z <- Z
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

    if(n.boot>0) objF$seed<-employed.Random.seed
    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last) list.obj <- list.obj[[length(list.obj)]]
#    warning("'segmented.Arima' is at a preliminary stage. Estimates are OK, but the '*.segmented' methods are not expected to work",
#      call.=FALSE)
    return(list.obj)
    } #end function
