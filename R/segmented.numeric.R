`segmented.numeric` <-
function(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(), model = TRUE, keep.class=FALSE, 
         adjX=FALSE, weights=NULL,  ...) { #sparse=FALSE,
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
  
  if(!(is.vector(obj) || is.ts(obj))) stop(" 'obj' should be a numerical/ts vector ")
  #if(!is.vector(obj)) stop(" 'obj' should be a numerical vector ")
  
  y <- obj
  n <- length(y)
  
  #browser()
  
  if(missing(seg.Z)) {
    if(is.ts(obj)){
      Tsp<-tsp(obj)
      x<-seq(Tsp[1], Tsp[2], length=length(y) )
      min.x<- min(x)
      name.Z <- "Time"
      if(is.null(adjX)) {
        adjX<- if(min.x>=1000) TRUE else FALSE
      } 
      if(adjX) x<- x - min.x
    } else {
      x<-1:n/n
      name.Z <- "index"
      adjX<-FALSE
    }
  } else {
    x<-eval(parse(text=all.vars(seg.Z)))
    name.Z <- all.vars(seg.Z)
    adjX= FALSE
  }
  
  #browser()
  
  if(!missing(seg.Z) && length(all.vars(seg.Z))>1) stop(" multiple seg.Z not allowed here: use 'segmented.(g)lm or segreg")
  Fo0<-as.formula(paste(deparse(substitute(obj)), " ~ ", name.Z, sep=""))
  #Fo0<-as.formula(paste(deparse(substitute(obj, env = parent.frame())), " ~ ", name.Z, sep="")) #qui mantiene il nome 'obj'
  
  #browser()
  y.only.vector <- TRUE
  alpha<-control$alpha
  if(is.null(alpha)) alpha<- max(.05, 1/length(y)) 
  if(length(alpha)==1) alpha<-c(alpha, 1-alpha)
  #browser()
  if(missing(psi)){
    if(missing(npsi)) npsi<-1 #stop(" psi or npsi have to be provided ")
    #psi <- seq(min(x), max(x), l=npsi+2)[-c(1, npsi+2)] #psi[[i]]<-(min(Z[[i]])+ diff(range(Z[[i]]))*(1:K)/(K+1))
    qx <- quantile(x, probs=c(alpha, 1-alpha), names = FALSE)
    psi <- seq(qx[1], qx[2], l=npsi+2)[-c(1, npsi+2)]
    
  } else {
    npsi<-length(psi)
  }
  a<- npsi
  initial.psi<-psi
  Z <- matrix(x, ncol=npsi,  nrow=n, byrow = FALSE)
  XREG <- cbind(1,x)
  PSI<-matrix(psi, ncol=a, nrow=n, byrow = TRUE)
  nomiU<-paste("U", 1:a, ".", name.Z,sep="")
  nomiV<-paste("V", 1:a, ".", name.Z,sep="")
  colnames(Z)<-nomiZ<-rep(name.Z, a)
  id.psi.group <- rep(1:length(a), times = a)
  orig.call<-NULL
  dev0<- n*var(y) #sum(mylm(x.lin, y, ww)$residuals^2*ww)
  rangeZ <- apply(Z, 2, range)
  if(is.null(weights)) {
    id.weights <- FALSE
    weights<-rep(1,n)
  } else {
    id.weights <- TRUE
  }

    fc<- min(max(abs(control$fc),.8),1)       
    min.step<-control$min.step
    
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
    
    break.boot=control$break.boot
    n.boot<-control$n.boot
    size.boot<-control$size.boot
    gap<-control$gap
    random<-control$random
    pow<-control$pow
    conv.psi<-control$conv.psi
    visual <- control$visual
    visualBoot<-FALSE
    if(visual && n.boot>0) {visual<-FALSE; visualBoot<-TRUE}
    
#     if(n.boot>0){
#         if(!is.null(control$seed)) {
#             set.seed(control$seed)
#             employed.Random.seed<-control$seed
#               } else {
#             employed.Random.seed<-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
#             set.seed(employed.Random.seed)
#               }
#         if(visual) {visual<-FALSE; visualBoot<-TRUE}# warning("`display' set to FALSE with bootstrap restart", call.=FALSE)}
# #        if(!stop.if.error) stop("Bootstrap restart only with a fixed number of breakpoints")
#      }
    last <- control$last
    K<-control$K
    h<-control$h
    #============================================
    #  ATTENZIONE devi costruire il mf con i pesi e la orig call?
    #============================================
    invXtX=NULL
    Xty<-NULL
    nomiOK<-nomiU
    
    opz<-list(toll=toll,h=h, stop.if.error=stop.if.error, dev0=dev0, visual=visual, it.max=it.max,
        nomiOK=nomiOK, id.psi.group=id.psi.group, gap=gap, visualBoot=visualBoot, pow=pow, digits=digits,invXtX=invXtX, Xty=Xty, 
        conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step, fc=fc, id.weights=id.weights, seed=control$seed, min.n=control$min.n)
    if(n.boot<=0){
    #obj<- if(sparse) seg.num.spar.fit(y, XREG, Z, PSI, weights,  opz) else seg.num.fit(y, XREG, Z, PSI, weights,  opz)
      obj<- seg.num.fit(y, XREG, Z, PSI, weights,  opz)
    } else {
    obj<-seg.num.fit.boot(y, XREG, Z, PSI, weights, opz, n.boot=n.boot, size.boot=size.boot, random=random, 
                          break.boot=break.boot) #, sparse=sparse) #jt, nonParam
    seed <- obj$seed
      }
    if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(y)
        }
    
    #browser()
    
    #if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
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
    idU <- obj$idU
    idV <- obj$idV
    obj<-obj$obj
    k<-length(psi)
    beta.c<-coef(obj)[idU]# [paste("U", 1:ncol(U), sep = "")]
    #psi.values[[length(psi.values) + 1]] <- psi #non c'e' bisogno!
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))

    #se usi una procedura automatica devi cambiare ripetizioni, nomiU e nomiV, e quindi:
    length.psi<-tapply(as.numeric(as.character(names(psi))), as.numeric(as.character(names(psi))), length)
    
    #forma.nomiU <-function(xx,yy)paste("U",1:xx, ".", yy, sep="")
    #forma.nomiVxb <-function(xx,yy)paste("psi",1:xx, ".", yy, sep="")
    #nomiU   <- unlist(mapply(forma.nomiU, length.psi, name.Z)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
    #nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, name.Z))
    #nomiU   <- unlist(mapply(forma.nomiU, length.psi, nomiFINALI)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
    #nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, nomiFINALI))
    nomiVxb <- sub("U","psi", nomiU)

    
    #########========================= SE PSI FIXED
    psi.list<-vector("list", length=length(unique(nomiZ)))
    names(psi.list)<-unique(nomiZ)
    #names(psi)<-nomiZ #se e' una procedure automatica nomiZ puo essere piu lungo dei breakpoints "rimasti" 
    names(psi)<-rep(nomiFINALI, length.psi)
    for(i in names(psi.list)){
      psi.list[[i]]<-psi[names(psi)==i]
    }

    #browser()
    
    #mf<- model.frame(update.formula(Fo0, .~ x))
    mf<- data.frame(y,x)
    names(mf)<-all.vars(Fo0)
    for(i in 1:ncol(U)) {
      #mfExt[nomiU[i]]<-
      mf[nomiU[i]]<-U[,i]
      #mfExt[nomiVxb[i]]<-
      mf[nomiVxb[i]]<-Vxb[,i]
    }
    nnomi <- c(nomiU, nomiVxb)
    #Fo <- update.formula(Fo0, as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
    #se c'e' un "y[-1]", la seguente linea modifica il nome in "y"..
    Fo <- update.formula(Fo0, as.formula(paste(paste(all.vars(Fo0)[1]),"~.+", paste(nnomi, collapse = "+"))))
    #mf <-  eval(mf, parent.frame()) #forse NON serve.. mf c'e'..
    objF <-lm(Fo, data=mf, weights=weights)
    #browser()
    
    isNAcoef<-any(is.na(objF$coefficients))
    if(isNAcoef){
      if(stop.if.error) {
       cat("breakpoint estimate(s):", as.vector(psi),"\n")
       stop("at least one coef is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)", 
         call. = FALSE)
          } else {
        warning("some estimate is NA: too many breakpoints? 'var(hat.psi)' cannot be computed \n ..returning a 'lm' model", call. = FALSE)
        Fo <- update.formula(Fo0, as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
        #objF <- update(obj0, formula = Fo,  evaluate=TRUE, data = mf)
        objF <-lm(Fo, weights=weights, data=mf)
        names(psi)<-nomiVxb
        objF$psi<-psi
        return(objF)      
        }
    }

    #browser()
    
    
    objF$coefficients[names(objF$coefficients)] <- obj$coefficients #sostituisce tutti i coeff 
    objF$residuals<- as.numeric(obj$residuals)
    objF$fitted.values<- y- as.numeric(obj$residuals) #as.numeric(obj$fitted.values) #y-obj$residuals 
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
    #browser()
    
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
    objF$orig.call<- update(objF, Fo0, evaluate=FALSE)
    objF$indexU<-build.all.psi(psi.list, fixed.psi)
    objF$psi[,"Initial"]<-NA
    if(model)  objF$model <- mf #objF$mframe <- data.frame(as.list(KK))
    if(n.boot>0) objF$seed<-seed
    class(objF) <- c("segmented", "lm")
    return(objF)
}
