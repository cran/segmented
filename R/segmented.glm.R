`segmented.glm` <-
#objF$id.group???
function(obj, seg.Z, psi=stop("provide psi"), control = seg.control(), model = TRUE, ...) {
    n.Seg<-1
    if(is.list(psi)){
      if(length(all.vars(seg.Z))!=length(psi)) stop("A wrong number of terms in `seg.Z' or `psi'")
      if(any(is.na(match(all.vars(seg.Z),names(psi), nomatch = NA)))) stop("Variables in `seg.Z' and `psi' do not match")
      n.Seg <- length(psi)
      }
    if(length(all.vars(seg.Z))!=n.Seg) stop("A wrong number of terms in `seg.Z' or `psi'")
    maxit.glm <- control$maxit.glm
    it.max <- old.it.max<- control$it.max
    toll <- control$toll
    visual <- control$visual
    stop.if.error<-control$stop.if.error
    n.boot<-control$n.boot
    size.boot<-control$size.boot
    gap<-control$gap
    visualBoot<-FALSE
    if(n.boot>0){
        if(visual) {visual<-FALSE; visualBoot<-TRUE}#warning("`display' set to FALSE with bootstrap restart", call.=FALSE)}
        if(!stop.if.error) stop("Bootstrap restart only with a fixed number of breakpoints")
     }
    last <- control$last
    K<-control$K
    h<-min(abs(control$h),1)
    if(h<1) it.max<-it.max+round(it.max/2)
#    if(!stop.if.error) objInitial<-obj        
    #-------------------------------
#    #una migliore soluzione.........
#    objframe <- update(obj, model = TRUE, x = TRUE, y = TRUE)
#    y <- objframe$y
#    a <- model.matrix(seg.Z, data = eval(obj$call$data))
#    a <- subset(a, select = colnames(a)[-1])
    orig.call<-Call<-mf<-obj$call
    orig.call$formula<-mf$formula<-formula(obj) #per consentire lm(y~.)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if(class(mf$formula)=="name" && !"~"%in%paste(mf$formula)) mf$formula<-eval(mf$formula)
    #orig.call$formula<-update.formula(orig.call$formula, paste("~.-",all.vars(seg.Z))) #utile per plotting
    
    nomeRispo<-strsplit(paste(formula(obj))[2],"/")[[1]] #eventuali doppi nomi (tipo "y/n" per GLM binom)
    #la linea sotto aggiunge nel mf anche la variabile offs..
    if(length(all.vars(formula(obj)))>1){
      id.rispo<-1
      if(length(nomeRispo)>=2) id.rispo<-1:2      
      #questo serve quando formula(obj) ha solo l'intercept
      agg<-if(length(all.vars(formula(obj))[-id.rispo])==0) "" else "+"
      mf$formula<-update.formula(mf$formula,paste(paste(seg.Z,collapse=".+"),agg,paste(all.vars(formula(obj))[-id.rispo],collapse="+")))
    } else {
      mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
    }
    mf <- eval(mf, parent.frame())
    #id.offs<-pmatch("offset",names(mf)) #questa identifica il nome offset(..). ELiminarlo dal dataframe? non conviene
    #       altrimenti nel model.frame non risulta l'offset
    weights <- as.vector(model.weights(mf))
    offs <- as.vector(model.offset(mf))
    
    if(!is.null(Call$weights)){ #"(weights)"%in%names(mf)
      names(mf)[which(names(mf)=="(weights)")]<-as.character(Call$weights)
      # mf["(weights)"]<-weights
      }
    
    mt <- attr(mf, "terms")
    interc<-attr(mt,"intercept")
    y <- model.response(mf, "any")
    XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    namesXREG0<-colnames(XREG)
    nameLeftSlopeZero<-setdiff(all.vars(seg.Z), all.vars(formula(obj)))
    namesXREG0<-setdiff(namesXREG0, nameLeftSlopeZero)
    
    #nomeRispo<-strsplit(paste(formula(obj))[2],"/")[[1]] #portato sopra
    if(length(nomeRispo)>=2) mf[nomeRispo[1]]<-weights*y
    
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
    #XREG<-model.matrix(obj0) non va bene perché non elimina gli eventuali mancanti in seg.Z..
    #Due soluzioni
    #XREG<-XREG[,colnames(model.matrix(obj)),drop=FALSE]
    #XREG<-XREG[,match(c("(Intercept)",all.vars(formula(obj))[-1]),colnames(XREG),nomatch =0),drop=FALSE]
    XREG <- XREG[, match(c("(Intercept)", namesXREG0),colnames(XREG), nomatch = 0), drop = FALSE]
    n <- nrow(XREG)
    #Z <- list(); for (i in colnames(XREGseg)) Z[[length(Z) + 1]] <- XREGseg[, i]
    Z<-lapply(apply(XREGseg,2,list),unlist) #prende anche i nomi!
    name.Z <- names(Z) <- colnames(XREGseg)
    if(length(Z)==1 && is.vector(psi) && (is.numeric(psi)||is.na(psi))){
        psi <- list(as.numeric(psi))
        names(psi)<-name.Z
        }
    if (!is.list(Z) || !is.list(psi) || is.null(names(Z)) || is.null(names(psi))) stop("Z and psi have to be *named* list")
    id.nomiZpsi <- match(names(Z), names(psi))
    if ((length(Z)!=length(psi)) || any(is.na(id.nomiZpsi))) 
        stop("Length or names of Z and psi do not match")
    #dd <- match(names(Z), names(psi))
    nome <- names(psi)[id.nomiZpsi]
    psi <- psi[nome]
    for(i in 1:length(psi)) {
        if(any(is.na(psi[[i]]))) psi[[i]]<-quantile(Z[[i]], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)
        }
    
    a <- sapply(psi, length)#b <- rep(1:length(a), times = a)
    id.psi.group <- rep(1:length(a), times = a) #identificativo di apparteneza alla variabile
    #Znew <- list()
    #for (i in 1:length(psi)) Znew[[length(Znew) + 1]] <- rep(Z[i], a[i])
    #Z <- matrix(unlist(Znew), nrow = n)
    Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n)
    psi <- unlist(psi)
    psi<-unlist(tapply(psi,id.psi.group,sort))
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
    colnames(Z) <- nomiZ <- rep(nome, times = a)
    ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))], function(xxx) {1:xxx})))
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
        U <- pmax((Z - PSI), 0)
        colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
        nomiU <- paste("U", colnames(U), sep = "")
        #for (i in 1:ncol(U)) assign(nomiU[i], U[, i], envir = KK)
        #è necessario il for? puoi usare colnames(U)<-nomiU;mf[nomiU]<-U
        for(i in 1:ncol(U)) mf[nomiU[i]]<-U[,i]
        Fo <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
        #obj <- update(obj, formula = Fo, data = KK)
        obj <- update(obj, formula = Fo, data = mf, evaluate=FALSE)
        if(!is.null(obj[["subset"]])) obj[["subset"]]<-NULL
        obj<-eval(obj, envir=mf)
        if (model) obj$model <-mf  #obj$model <- data.frame(as.list(KK))
        names(psi)<-paste(paste("psi", ripetizioni, sep = ""), nomiZ, sep=".")
        obj$psi <- psi
        return(obj)
    }
    #XREG <- model.matrix(obj) creata sopra         
    #o <- model.offset(objframe)
    #w <- model.weights(objframe)
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(offs)) offs <- rep(0, n)
    fam <- family(obj)
    initial <- psi
    obj0 <- obj
    dev0<-obj$dev
    list.obj <- list(obj)
#    psi.values <- NULL
    nomiOK<-nomiU
    opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,dev0=dev0,visual=visual,it.max=it.max,nomiOK=nomiOK,
        fam=fam, eta0=obj$linear.predictors, maxit.glm=maxit.glm, id.psi.group=id.psi.group,gap=gap,visualBoot=visualBoot)

    if(n.boot<=0){
      obj<-seg.glm.fit(y,XREG,Z,PSI,weights,offs,opz)
    } else {
      obj<-seg.glm.fit.boot(y, XREG, Z, PSI, weights, offs, opz, n.boot=n.boot, size.boot=size.boot)
      }
    if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(obj0)
        }
    nomiOK<-obj$nomiOK
    it<-obj$it
    psi<-obj$psi
    k<-length(psi)
    psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
    U<-obj$U
    V<-obj$V
    rangeZ<-obj$rangeZ 
    obj<-obj$obj
    beta.c<-if(k == 1) coef(obj)["U"] else coef(obj)[paste("U", 1:ncol(U), sep = "")]
    psi.values[[length(psi.values) + 1]] <- psi
    id.warn <- FALSE
    if (n.boot<=0 && it > it.max) {
        warning("max number of iterations attained", call. = FALSE)
        id.warn <- TRUE
    }
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
    colnames(U) <- colnames(Vxb) <-sapply(strsplit(nomiOK,"U"),function(x)x[2])
    #colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
    #colnames(Vxb) <- paste(ripetizioni, nomiZ, sep = ".")
    nomiU <- paste("U", colnames(U), sep = "")
    nomiVxb <- paste("psi", colnames(Vxb), sep = "")
    for(i in 1:ncol(U)) {
        mf[nomiU[i]]<-U[,i]
        mf[nomiVxb[i]]<-Vxb[,i]
        }
#    for (i in 1:ncol(U)) {
#        assign(nomiU[i], U[, i], envir = KK)
#        assign(nomiVxb[i], Vxb[, i], envir = KK)
#    }
    nnomi <- c(nomiU, nomiVxb)
    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", 
        paste(nnomi, collapse = "+"))))
    if(is.matrix(y)&& (fam$family=="binomial" || fam$family=="quasibinomial")){
              mf<-cbind(mf[[1]], mf[,-1])
    }
    objF <- update(obj0, formula = Fo, data = mf, evaluate=FALSE)
    if(!is.null(objF[["subset"]])) objF[["subset"]]<-NULL
    objF<-eval(objF, envir=mf)
#c'è un problema..controlla obj (ha due "(Intercepts)"
    if(!gap){
        names.coef<-names(objF$coefficients)
        objF$coefficients<- if(sum("(Intercept)"==names(obj$coef))==2) obj$coefficients[-2] else obj$coefficients
        names(objF$coefficients)<-names.coef
        objF$fitted.values<-obj$fitted.values
        objF$residuals<-obj$residuals
        }
    if(any(is.na(objF$coefficients))){
     stop("some estimate is NA: premature stopping with a large number of breakpoints?",
      call. = FALSE)
    }
    Cov <- vcov(objF)
    id <- match(nomiVxb, names(coef(objF)))
    #cat(id,"\n")
    #return(objF)
    vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
    if(length(initial)!=length(psi)) initial<-rep(NA,length(psi)) 
    psi <- cbind(initial, psi, sqrt(vv))
    rownames(psi) <- colnames(Cov)[id]
    colnames(psi) <- c("Initial", "Est.", "St.Err")
    objF$rangeZ <- rangeZ
    objF$psi.history <- psi.values
    objF$psi <- psi
    objF$it <- (it - 1)
    objF$epsilon <- obj$epsilon
    objF$call <- match.call()
    objF$nameUV <- list(U = nomiU, V = rownames(psi), Z = name.Z)
    objF$id.group <- if(length(name.Z)<=1) -rowSums(as.matrix(V))
    objF$id.warn <- id.warn
    objF$orig.call<-orig.call
    if (model)  objF$model <- mf #objF$mframe <- data.frame(as.list(KK))
    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last) 
        list.obj <- list.obj[[length(list.obj)]]
    return(list.obj)
}
