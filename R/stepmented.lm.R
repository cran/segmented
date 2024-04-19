stepmented.lm <- function(obj, seg.Z, psi, npsi, fixed.psi=NULL, control=seg.control(), 
                          keep.class=FALSE, var.psi=FALSE, ...) {
  # ---------
  mylm<-function(x,y,w=1,offs=0){
    x1<-x*sqrt(w)
    y<-y-offs
    y1<-y*sqrt(w)
    XtX <- crossprod(x1)
    b<-drop(solve(XtX,crossprod(x1,y1)))
    fit<-drop(tcrossprod(x,t(b)))
    r<-y-fit
    o<-list(coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b), invXtX=solve(XtX), w=w)
    o
  }
  #-----------
  toMatrix<-function(.x, ki){
    # ripete ogni .x[,j] ki[j] volte
    if(ncol(.x)!=length(ki)) stop("It should be ncol(.x)==length(ki)")
    if(all(ki==1)) return(.x)
    M<-vector("list", length=length(ki))
    for(j in 1:length(ki)) M[[j]]<-replicate(ki[[j]], cbind(.x[,j]), simplify=TRUE)
    do.call(cbind, M)
  }
  #-----------
  agg<- 1-control$fc
  it.max<- control$it.max 
  tol<-  control$toll
  display<- control$visual 
  digits <- control$digits 
  min.step <- control$min.step
  #conv.psi <- control$conv.psi 
  alpha <- control$alpha
  fix.npsi <- control$fix.npsi 
  n.boot <- control$n.boot 
  break.boot<- control$break.boot + 2 
  seed<- control$seed
  fix.npsi<-control$fix.npsi
  h<-control$h
  #-----------
  #browser()
  
  #if(!(inherits(obj,"lm") || is.vector(obj) || is.ts(obj))) stop("obj should be a 'lm' fit, a 'vector' or 'ts' object")
  if(!inherits(obj,"lm")) stop("obj should be a 'lm' fit")
  
    y.only.vector <- FALSE
    Fo0 <- formula(obj) 
    if(missing(seg.Z)) {
      #if(length(all.vars(formula(obj)))==1) 
      seg.Z<- as.formula(paste("~", "id"))
      assign("id",1:length(obj$residuals),parent.frame())
      #id<-1:length(obj$residuals)
      #    if(length(all.vars(formula(obj)))==2) seg.Z<- as.formula(paste("~", all.vars(formula(obj))[2])) else stop("please specify 'seg.Z'")
    }
    n.Seg<-length(all.vars(seg.Z))
    id.npsi<-FALSE
    
    if("V" %in% sub("V[1-9]*[0-9]","V", c(all.vars(seg.Z), all.vars(formula(obj) )[-1]))) stop("variable names 'V', 'V1', .. are not allowed")
    if("U" %in% sub("U[1-9]*[0-9]","U", c(all.vars(seg.Z), all.vars(formula(obj) )[-1]))) stop("variable names 'U', 'U1', .. are not allowed")
    if(any(c("$","[") %in% all.names(seg.Z))) stop(" '$' or '[' not allowed in 'seg.Z' ")
    
    if(missing(psi)){ 
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
    n.psi<- length(unlist(psi))
    
    #browser()
    
    #if(missing(x)) x<-1:n
    #if(missing(psi)) {
    #  if(missing(npsi)) npsi<-1
    #  psi<-min(x)+cumsum(rep(diff(range(x))/(npsi+1),npsi))
    #}
    ##=========================================================================
    #--- preso da segmented.lm
    orig.call<-Call<-mf<-obj$call
    orig.call$formula<- mf$formula<-formula(obj) #per consentire lm(y~.)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if(class(mf$formula)[1]=="name" && !"~"%in%paste(mf$formula)) mf$formula<-eval(mf$formula)
    mfExt<- mf
    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
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
    nomiTUTTI<-all.vars(mfExt$formula) #comprende anche altri nomi (ad es., threshold) "variabili"
    nomiNO<-NULL 
    for(i in nomiTUTTI){
      r<-try(eval(parse(text=i), parent.frame()), silent=TRUE)
      if(class(r)[1]!="try-error" && length(r)==1 && !is.function(r) && !i%in%names(mf)) nomiNO[[length(nomiNO)+1]]<-i
    }
    if(!is.null(nomiNO)) mfExt$formula<-update.formula(mfExt$formula,paste(".~.-", paste( nomiNO, collapse="-"), sep=""))
    mfExt<-eval(mfExt, parent.frame())
    
    #mf <- mfExt
    #browser()
    if(nrow(mf)!=nrow(mfExt)) stop("missing values in any stepmented covariate?")
    
    
    ww <- as.vector(model.weights(mf))
    offs <- as.vector(model.offset(mf))
    if (is.null(ww)) ww <- rep(1, n)
    if (is.null(offs)) offs <- rep(0, n)
    mt <- attr(mf, "terms")
    interc<-attr(mt,"intercept")
    y <- model.response(mf, "any")
    XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf, obj$contrasts)
    namesXREG0<-colnames(XREG)
    nameLeftSlopeZero<-setdiff(all.vars(seg.Z), names(coef(obj))) #in questo modo riconosce che sin(x*pi) NON e' x, ad esempio.
    namesXREG0<-setdiff(namesXREG0, nameLeftSlopeZero)
    id.duplic<-match(all.vars(formula(obj)),all.vars(seg.Z),nomatch=0)>0
    if(any(id.duplic)) {
      new.mf<-mf[,all.vars(formula(obj))[id.duplic],drop=FALSE]
      new.XREGseg<-data.matrix(new.mf)
      XREG<-cbind(XREG,new.XREGseg)
    }
    
    id.n.Seg<-(ncol(XREG)-n.Seg+1):ncol(XREG)
    XREGseg<-XREG[,id.n.Seg,drop=FALSE]
    XREG <- XREG[, match(c("(Intercept)", namesXREG0),colnames(XREG), nomatch = 0), drop = FALSE]
    XREG<-XREG[,unique(colnames(XREG)), drop=FALSE]
    n <- nrow(XREG)
    
    #browser()

    Z<-lapply(apply(XREGseg,2,list),unlist) #prende anche i nomi!
    name.Z <- names(Z) <- colnames(XREGseg)
    
    if(length(Z)==1 && is.vector(psi) && (is.numeric(psi)||is.na(psi))){
      psi <- list(as.numeric(psi))
      names(psi)<-name.Z
    }
    if (!is.list(Z) || !is.list(psi) || is.null(names(Z)) || is.null(names(psi))) stop("'psi' or 'npsi' have to be *named* when there are multiple stepmented variables")
    id.nomiZpsi <- match(names(Z), names(psi))
    if ((length(Z)!=length(psi)) || any(is.na(id.nomiZpsi))) stop("Length or names of 'seg.Z' and 'psi' do not match")
    nome <- names(psi)[id.nomiZpsi]
    psi <- psi[nome]
    if(id.npsi){
      for(i in 1:length(psi)) {
        K<-length(psi[[i]])
        if(any(is.na(psi[[i]]))) psi[[i]]<-(min(Z[[i]])+ diff(range(Z[[i]]))*(1:K)/(K+1))
      }
    } else {
      for(i in 1:length(psi)) {
        if(any(is.na(psi[[i]]))) psi[[i]]<- (min(Z[[i]])+ diff(range(Z[[i]]))*(1:K)/(K+1))
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
    a <- sapply(psi, length) #n. di psi per ogni covariate
    #per evitare che durante il processo iterativo i psi non siano ordinati
    id.psi.group <- rep(1:length(a), times = a) #identificativo di apparteneza alla variabile
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
    
    #browser()
    
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ, sep = ".")
    nomiV <- paste("V", ripetizioni, sep = "")
    nomiV <- paste(nomiV, nomiZ, sep = ".")
    initial <- psi
    obj0 <- obj
    dev0 <-sum(ww*obj$residuals^2)
    list.obj <- list(obj)
    nomiOK<-nomiU
  
  #  invXtX<-if(!is.null(obj$qr)) chol2inv(qr.R(obj$qr)) else NULL #(XtX)^{-1}
  #  Xty<-crossprod(XREG,y)
  #  opz<-list(toll=toll,h=h, stop.if.error=stop.if.error, dev0=dev0, visual=visual, it.max=it.max,
  #            nomiOK=nomiOK, id.psi.group=id.psi.group, gap=gap, visualBoot=visualBoot, pow=pow, digits=digits,invXtX=invXtX, Xty=Xty, 
  #            conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step, fc=fc)
  x.lin <-XREG
  rangeZ <- apply(Z, 2, range)
  
  #browser()
  
  plin<-ncol(x.lin)
  #if(!is.list(psi)) psi<-list(psi)
  #P <- length(psi) #n. variabili con breakpoints
  #npsii <- sapply(psi, length) #n di psi for each covariate
  P<-n.Seg
  npsii<-a
  npsi<- sum(npsii)
  #Xtrue<-Z
  #psi0 <- unlist(psi)
  #PSI<- matrix(psi0, n, npsi, byrow=TRUE)
  #if(ncol(x)!=P) stop("errore")
  #Xtrue<-toMatrix(x, npsii)
  
  #browser()
  
  if(it.max == 0) {
    U <- (Z>PSI)
    colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
    nomiU <- paste("U", colnames(U), sep = "")
    #for (i in 1:ncol(U)) assign(nomiU[i], U[, i], envir = KK)
    for(i in 1:ncol(U)) mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
    Fo <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
    obj <- update(obj, formula = Fo, evaluate=FALSE, data=mfExt) #data = mf, 
    if(!is.null(obj[["subset"]])) obj[["subset"]]<-NULL
    obj<-eval(obj, envir=mfExt)
    #if (model) obj$model <-mf  #obj$model <- data.frame(as.list(KK))
    
    psi <- cbind(psi, psi, 0)
    rownames(psi) <- paste(paste("psi", ripetizioni, sep = ""), nomiZ, sep=".")
    colnames(psi) <- c("Initial", "Est.", "St.Err")
    
    obj$psi <- psi
    return(obj)
  }
  
  
  c1 <- apply((Z <= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo <)
  c2 <- apply((Z >= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo >)
  if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2)) ) stop("starting psi out of the admissible range")
  if(is.null(alpha)) alpha<- max(.05, 1/length(y))
  if(length(alpha)==1) alpha<-c(alpha, 1-alpha)
  
  #browser()
  
  opz<-list(toll=tol, dev0=dev0, display=display, it.max=it.max, agg=agg, digits=digits, rangeZ=rangeZ, usestepreg=FALSE,
            id.psi.group=id.psi.group, h=h,
            #nomiOK=nomiOK, , visualBoot=visualBoot, invXtX=invXtX, Xty=Xty, conv.psi=conv.psi, 
            alpha=alpha, fix.npsi=fix.npsi, min.step=min.step, npsii=npsii, 
            seed=control$seed, fit.psi0=control$fit.psi0)
  
  # #################################################################################
  # #### jump.fit(y, XREG=x.lin, Z=Xtrue, PSI, w=ww, offs, opz, return.all.sol=FALSE)
  # #################################################################################
  if(n.boot<=0){
    obj<- step.lm.fit(y, x.lin, Z, PSI, ww, offs, opz, return.all.sol=FALSE)
  } else {
    #browser()
    #if("seed" %in% names(control)) set.seed(control$seed)
    obj<-step.lm.fit.boot(y, x.lin, Z, PSI, ww, offs, opz, n.boot, break.boot=break.boot) 
    seed<- obj$seed
  }
  # if(!is.list(obj)){
  #   warning("No breakpoint estimated", call. = FALSE)
  #   return(obj0)
  # }

  #browser()
  id.warn<-obj$id.warn
  it<-obj$it
  psi<-obj$psi
  psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
  #i beta.c corripondono ai psi NON ordinati!!!
  #
  ##Nelle funzioni step i beta.c NON servono. Righe sotto commentate l'8/3/24
  #beta.c<- obj$beta.c
  #beta.c<-unlist(tapply(psi, id.psi.group, function(.x)beta.c[order(.x)]))
  #unlist(lapply(unique(id.psi.group), function(.x) beta.c[id.psi.group==.x][order(psi[id.psi.group==.x])]))
  psi<-unlist(tapply(psi, id.psi.group, sort)) 
  Z0<-apply(Z,2,sort)
  psi.rounded<-sapply(1:npsi, function(j) Z0[sum(Z0[,j]<psi[j])+c(0,1),j])
  
  #browser()
  
  psi.mid<-apply(psi.rounded,2,mean)
  #QUALI prendere? psi, psi.mid o psi.rounded?
  PSI.mid<- matrix(psi, n, npsi, byrow = TRUE)
  
  #bisogna evitare che una qualche x_i sia uguale a psi, altrimenti la costruzione di V-> INF
  DEN <- abs(Z - PSI.mid)
  DEN <- apply(DEN, 2, function(.x) pmax(.x, sort(.x)[2]/2))  #pmax(.x, diff(range(.x))/1000)) 
  
  #xx=Xtrue - PSI.mid
  #ss=n^(-.8)
  #den <- -xx+2*xx*pnorm(xx/ss)+2*ss*dnorm(xx/ss)     #.05*log(cosh((x-.5)/.05)))
  
  V <- (1/(2 * DEN))
  colnames(V)<-nomiV
  U <- (Z * V + 1/2)
  colnames(U)<-nomiU
  Vxb <- -V  #* rep(-beta.c, each = nrow(V))
  nomiVxb <- gsub("V", "psi", nomiV)
  nnomi <- c(nomiU, nomiVxb)
  #browser()
  
  for(i in 1:ncol(U)) {
      mfExt[nomiU[i]]<-mf[nomiU[i]] <- U[,i]
      mfExt[nomiVxb[i]]<-mf[nomiVxb[i]] <- Vxb[,i]
  }
  Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
  objF <- update(obj0, formula = Fo,  evaluate=FALSE, data = mfExt)
  #eliminiamo subset, perche' se e' del tipo subset=x>min(x) allora continuerebbe a togliere 1 osservazione 
  if(!is.null(objF[["subset"]])) objF[["subset"]]<-NULL
  objF<-eval(objF, envir=mfExt)
  objF$offset<- obj0$offset
  objW<-objF
  
  #browser()
  
  #se1=predict.lm(objF, se.fit=TRUE)
  #ff<-1.934+1.61*(x>.605)
  #matplot(x, cbind(ff, ff-2*se$se.fit, ff+2*se$se.fit), type="l")
  
  #controllo se qualche coeff e' NA..
  isNAcoef<-any(is.na(objF$coefficients))
 
 #browser()
  if (isNAcoef) {
    nameNA.psi <- names(objF$coefficients)[which(is.na(objF$coefficients))]
    nameNA.U <- gsub("psi", "U", nameNA.psi)
    if (fix.npsi) {
        cat("breakpoint estimate(s):", as.vector(psi), "\n")
        stop("coef ", nameNA.psi, " is NA: breakpoint(s) at the boundary or too close together", call. = FALSE)
    } else {
        warning("some estimate is NA (too many breakpoints?): removing ", length(nameNA.psi), " jump-point(s)", call. = FALSE)
        Fo <- update(Fo, paste(".~ .-", nameNA.U, "-", nameNA.psi))
        objF <- update(obj0, formula = Fo, evaluate = TRUE, data = mfExt)
        if (!is.null(objF[["subset"]])) objF[["subset"]] <- NULL
        #objF$offset <- obj0$offset
        idNA.psi <- match(nameNA.psi, nomiVxb)
        nomiVxb <- setdiff(nomiVxb, nameNA.psi)
        nomiU <- setdiff(nomiU, nameNA.U)
        Z <- Z[, -idNA.psi, drop = FALSE]
        PSI.mid<- PSI.mid[, -idNA.psi, drop = FALSE]
        id.psi.group <- id.psi.group[-idNA.psi]
        psi <- psi[-idNA.psi]
        psi.rounded <- psi.rounded[, -idNA.psi, drop = FALSE]
    }
  }
  
  #organizziamo i risultati da restituire per psi...
  colnames(psi.rounded)<-names(psi)<-nomiVxb
  rownames(psi.rounded)<-c("inf [","sup (")
  

  #browser()
  
  
  ris.psi<-matrix(NA,length(psi), 3)
  colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
  rownames(ris.psi) <- nomiVxb
  ris.psi[,2]<-psi
  #ris.psi[,3]<-sqrt(vv)
  
  a<-tapply(id.psi.group, id.psi.group, length)
  #NB "a" deve essere un vettore che si appatta con "initial.psi" per ottnetere "initial" sotto... Se una variabile alla fine risulta
  # senza breakpoint questo non avviene e ci sono problemi nella formazione di "initial". Allora costruisco a.ok
  a.ok<-NULL
  nomiFINALI<-unique(nomiZ)
  
  for(j in name.Z){
    if(j %in% nomiFINALI) { 
      a.ok[length(a.ok)+1]<-a[1]
      a<-a[-1]
    } else {
      a.ok[length(a.ok)+1]<-0
    } #ifelse(name.Z %in% nomiFINALI,1,0)
  }
  #initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi, a.ok, SIMPLIFY = TRUE))
  if(length(psi)!=length(initial.psi)){
    ris.psi[,1]<- NA
  } else {
    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi[nomiFINALI], a.ok[a.ok!=0], SIMPLIFY = TRUE))
    ris.psi[,1]<-initial #if(stop.if.error)  ris.psi[,1]<-initial 
  }
 
  objF$psi <- ris.psi
  objF$psi.rounded <- psi.rounded
  #objW<-objF
  #stima il modello "vero" (non-working)
  U<- (Z > PSI.mid)
  colnames(U)<-nomiU
  X <- cbind(x.lin, U)
  objF$objW<- objW
  objF$obj.ok<-mylm(X, y, w=ww, offs=offs) #coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b))
  objF$fitted.values<-objF$obj.ok$fitted.values
  objF$residuals<- objF$obj.ok$residuals
  objF$coefficients[1:length(objF$obj.ok$coefficients)] <- objF$obj.ok$coefficients
  objF$coefficients[nomiVxb] <-psi.rounded[1,]
  objF$nameUV <- list(U = drop(nomiU), V = nomiV, Z = name.Z) #Z = name.Z
  objF$rangeZ<-obj$rangeZ
  objF$Z<-Z[,unique(name.Z),drop=FALSE]
  objF$call <- match.call()
  objF$orig.call<-orig.call
  objF$psi.history <- psi.values
  objF$it <- it 
  objF$epsilon <- obj$epsilon
  objF$id.warn <- id.warn
  if(n.boot>0) objF$seed <- seed
  class(objF) <- c("stepmented", class(obj0))
  
  #Un effetto aggiuntivo..
  Z.in.obj<-intersect(all.vars(Fo0), all.vars(seg.Z))
  if(length(Z.in.obj)>0){
    tt<-terms(Fo0)#, specials=Z.in.obj)
    #id<-match(Z.in.obj, all.vars(Fo0))-1 #1 e' per la risposta..
    id<-match(Z.in.obj, intersect(all.vars(Fo0), names(mf)))-1
      
    nome<-attr(tt,"term.labels")[id]
    Fo.ok<-as.formula(paste("~0", nome, sep="+"))
    f.x<-matrix(NA, 150, ncol(objF$Z[,Z.in.obj,drop=FALSE])) #prima era nrow(objF$Z) invece che 100
    for(j in 1:length(Z.in.obj)){
      idPsi <- nomiVxb[endsWith(nomiVxb, paste(".", Z.in.obj[j], sep = ""))]
      #psi <- coef(objF)[idPsi]
      dd<-data.frame(seq(min(objF$Z[,Z.in.obj[j]]), max(objF$Z[,Z.in.obj[j]]), l=nrow(f.x)))
      names(dd)<- Z.in.obj[j]
      M<-model.matrix(Fo.ok, data=dd)
      f.x[,j]<-M%*% coef(objF)[colnames(M)]
    }
    colnames(f.x)<-Z.in.obj
    objF$f.x<-f.x
  }
  
  objF$psi<- objF$psi[,-1,drop=FALSE] #rimuovi la colonna Initial
  
  if(var.psi){
    Cov <- vcov.stepmented(objF, k=NULL)
    id <- match(nomiVxb, names(coef(objF)))
    vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
    objF$psi[,"St.Err"]<-sqrt(vv)
    objF$vcov<- Cov
  }
  #Cov[nomiVxb, ]<- Cov[, nomiVxb] <- 0

  # var.Tay<-function(est1,est2,v1,v2,v12){
  #   r<- est1/est2
  #   vv<-(v1+v2*r^2-2*r*v12)/est2^2
  #   vv}
  # 
  # varPsi<- rep(NA, length(nomiU))
  # for(j in 1:length(nomiU)){
  #   num<-objF$coefficients[nomiVxb[j]]
  #   den<-objF$coefficients[nomiU[j]]
  #   v.g <-Cov[nomiVxb[j],nomiVxb[j]]
  #   v.b<- Cov[nomiU[j],nomiU[j]]
  #   cov.g.b <- Cov[nomiVxb[j],nomiU[j]]
  #   #if(is.null(rho)) 
  #   rho<-mean(Xtrue[, nomiZ[j] ,drop=TRUE]<psi[[nomiVxb[j]]])
  #   #browser()
  #   rho<-  rho^(sqrt(1/n))
  #   cov.g.b<- rho*sqrt(v.g*v.b)
  #   varPsi[j]<-var.Tay(num, den, v.g, v.b, cov.g.b)
  # }
  # names(varPsi) <- nomiVxb

  return(objF)
}
