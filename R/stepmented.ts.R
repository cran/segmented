stepmented.ts <- function(obj, seg.Z, psi, npsi, fixed.psi=NULL, control=seg.control(), keep.class=FALSE, 
                          var.psi=FALSE, ..., pertV=0, centerX=FALSE, adjX=NULL) {
  #pertV come calcolare la variabile V=1/(2*abs(Xtrue-PSI)? i psi devono essere diversi dalle x_i 
  #   utilizzare i psi stimati che tipcamente sono diversi? (perV=0)
  #   oppure i psi.mid che sicuramente sono (o meglio dovrebbero essere)  tra due x_i...
  
  #NO OFFSET o PESI!!!!
  # ---------
  mylm<-function(x,y){
    XtX <- crossprod(x)
    b <- drop(solve(XtX,crossprod(x,y)))
    fit <- drop(tcrossprod(b,x))
    r<-y-fit
    o<-list(coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b), invXtX=solve(XtX), L0=sum(r^2))
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
  #browser()
  
  display<- control$visual 
  digits <- control$digits 
  #min.step <- control$min.step
  #conv.psi <- control$conv.psi 
  alpha <- control$alpha
  fix.npsi <- control$fix.npsi 
  n.boot <- control$n.boot 
  break.boot<- control$break.boot + 2
  seed<- control$seed
  fix.npsi<-control$fix.npsi
  h<-control$h
  #-----------
  #if(!(is.vector(obj) || is.ts(obj))) stop("obj should be a 'lm' fit, a 'vector' or 'ts' object")
  
  #if(is.vector(obj) || is.ts(obj)){
  #if(is.matrix(obj) && ncol(obj)>1) stop("if matrix 'obj' should have 1 column")
  #obj<-drop(obj)
  if(!missing(seg.Z) && length(all.vars(seg.Z))>1) stop(" multiple seg.Z allowed only with lm models")
  Fo0<-as.formula(paste(deparse(substitute(obj))," ~ 1", sep=""))
  y.only.vector<-TRUE
  y<- as.vector(obj)
  Tsp<-tsp(obj) 
  nn <- 1 + round((Tsp[2L] - Tsp[1L]) * Tsp[3L])
  #preso da getAnywhere("print.ts")

  if (length(y) != nn) warning(gettextf("series is corrupt: length %d with 'tsp' implying %d", length(y), nn), domain = NA, call. = FALSE)
  
  if(missing(seg.Z)) {
        x<-seq(Tsp[1], Tsp[2], length=length(y) )
        name.Z <- "Time"
  } else {
        x<-eval(parse(text=all.vars(seg.Z)))
        name.Z <- all.vars(seg.Z)
        adjX <- FALSE
  }
  min.x<- min(x)
  if(is.null(adjX)) {
    adjX<- if(min.x>=1000) TRUE else FALSE
  } 
  if(adjX) x<- x - min.x
  

  if(missing(psi)){
      if(missing(npsi)) npsi<-1 #stop(" psi or npsi have to be provided ")
      psi<- seq(min(x), max(x), l=npsi+2)[-c(1, npsi+2)] #psi[[i]]<-(min(Z[[i]])+ diff(range(Z[[i]]))*(1:K)/(K+1))
    } else {
      npsi<-length(psi)
  }
  initial.psi<-psi
  n<-length(y)
  a<- npsi
  n.Seg<-1
  Z<-matrix(x, ncol=a,  nrow=n, byrow = FALSE)
  XREG<-matrix(1, nrow=n, ncol=1)
  PSI<-matrix(psi, ncol=a, nrow=n, byrow = TRUE)
  #name.Z <- if(missing(seg.Z)) "id" else all.vars(seg.Z)
  nomiU<-paste("U", 1:a, ".", name.Z,sep="")
  nomiV<-paste("V", 1:a, ".", name.Z,sep="")
  colnames(Z)<-nomiZ<-rep(name.Z, a)
  id.psi.group <- rep(1:length(a), times = a)
  orig.call<-NULL
  ####################################################

  #  invXtX<-if(!is.null(obj$qr)) chol2inv(qr.R(obj$qr)) else NULL #(XtX)^{-1}
  #  Xty<-crossprod(XREG,y)
  #  opz<-list(toll=toll,h=h, stop.if.error=stop.if.error, dev0=dev0, visual=visual, it.max=it.max,
  #            nomiOK=nomiOK, id.psi.group=id.psi.group, gap=gap, visualBoot=visualBoot, pow=pow, digits=digits,invXtX=invXtX, Xty=Xty, 
  #            conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step, fc=fc)
  #x<- Z
  x.lin <-XREG
  #if(is.vector(x)) x<-as.matrix(x)
  #dev0<- if(!display) var(y)*n else sum(mylm(x.lin, y)$residuals^2)
  #non ci sono altre esplicative, per cui il modello nullo e' sempre con interc.
  dev0 <- var(y)*(n-1) #mylm(x.lin, y)$L0
  rangeZ <- apply(Z, 2, range)
  
  #browser()
  
  plin<-ncol(x.lin)
  #if(!is.list(psi)) psi<-list(psi)
  #P <- length(psi) #n. variabili con breakpoints
  #npsii <- sapply(psi, length) #n di psi for each covariate
  P<-n.Seg
  npsii<-a
  npsi<- sum(npsii)
  Xtrue<-Z
  #psi0 <- unlist(psi)
  #PSI<- matrix(psi0, n, npsi, byrow=TRUE)
  #if(ncol(x)!=P) stop("errore")
  #Xtrue<-toMatrix(x, npsii)
  
  #browser()
  
  if(it.max == 0) {
    mfExt<-data.frame(y, Z)
    names(mfExt) <- c(all.vars(Fo0), name.Z)
    ripetizioni<-unlist(tapply(nomiZ, nomiZ, function(.x)1:length(.x)))
    U <- (Xtrue>PSI)
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
  
  
  c1 <- apply((Xtrue <= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo <)
  c2 <- apply((Xtrue >= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo >)
  if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2)) ) stop("starting psi out of the admissible range")
  if(is.null(alpha)) alpha<- max(.05, 1/length(y))
  if(length(alpha)==1) alpha<-c(alpha, 1-alpha)
  
  opz<-list(toll=tol, dev0=dev0, display=display, it.max=it.max, agg=agg, digits=digits, rangeZ=rangeZ,
            id.psi.group=id.psi.group,h=h,
            #nomiOK=nomiOK,  visualBoot=visualBoot, invXtX=invXtX, Xty=Xty, conv.psi=conv.psi,min.step=min.step,
             alpha=alpha, fix.npsi=fix.npsi,  npsii=npsii,
            seed=control$seed)
  
  # #################################################################################
  # #### jump.fit(y, XREG=x.lin, Z=Xtrue, PSI, w=ww, offs, opz, return.all.sol=FALSE)
  # #################################################################################

  if(n.boot<=0){
    obj<- step.ts.fit(y, x.lin, Xtrue, PSI, opz, return.all.sol=FALSE)
  } else {
    #if("seed" %in% names(control)) set.seed(control$seed)
    obj<-step.ts.fit.boot(y, x.lin, Xtrue, PSI, opz, n.boot, break.boot=break.boot)
    seed <- obj$seed
  }
  # if(!is.list(obj)){
  #   warning("No breakpoint estimated", call. = FALSE)
  #   return(obj0)
  # }
  #chol2inv(qr.R(obj$obj$qr))
  id.warn<-obj$id.warn
  it<-obj$it
  psi<-obj$psi
  psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
  #i beta.c corripondono ai psi NON ordinati!!!
  beta.c<- obj$beta.c
  beta.c<-unlist(tapply(psi, id.psi.group, function(.x)beta.c[order(.x)]))
  #unlist(lapply(unique(id.psi.group), function(.x) beta.c[id.psi.group==.x][order(psi[id.psi.group==.x])]))
  psi<-unlist(tapply(psi, id.psi.group, sort)) 
  Z0<-apply(Z,2,sort)
  psi.rounded<-sapply(1:npsi, function(j) Z0[sum(Z0[,j]<psi[j])+c(0,1),j])
  psi.mid<-apply(psi.rounded,2,mean)
  #QUALI prendere? psi, psi.mid o psi.rounded?
  PSI.mid<- matrix(psi, n, npsi, byrow = TRUE)

  #bisogna evitare che una qualche x_i sia uguale a psi, altrimenti la costruzione di V-> INF
  DEN <- abs(Xtrue - PSI.mid)
  DEN <- apply(DEN, 2, function(.x) pmax(.x, sort(.x)[2]/2))  #pmax(.x, diff(range(.x))/1000)) 
  
  V <- (1/(2 * DEN))
  colnames(V)<-nomiV
  if(centerX){
    XtrueS <- scale(Xtrue, TRUE, scale=FALSE)
    meanX<-attr(XtrueS, "scaled:center")
    attr(XtrueS, "scaled:center")<-NULL
    U <- (XtrueS * V + 1/2)
  } else {
  U <- (Xtrue * V + 1/2)
  }
  colnames(U)<-nomiU

  if(pertV>0){
    #puoi usare o psi.mid o psi.rounded+eps.. Il secondo porta ad una cor ancora piu bassa della prima.. 0.89 vs 0.96
    if(pertV==1){
      PSI.mid <- matrix(psi.mid, n, npsi, byrow = TRUE)
      V <- (1/(2 * abs(Xtrue - PSI.mid)))
    } else {
      PSI.mid <- matrix(psi.rounded[1,], n, npsi, byrow = TRUE)
      V <- (1/(2 * abs(Xtrue - PSI.mid + .0001)))
    }
  }

  Vxb <- -V# * rep(-beta.c, each = nrow(V))
  nomiVxb <- gsub("V", "psi", nomiV)
  nnomi <- c(nomiU, nomiVxb)
  #XREG <- cbind(x.lin, Z, W)
  #obj <- lm.wfit(y = y, x = XREG, offset = offs, w=ww )
  # source("stepmented.lm.R")
  Fo <- update.formula(Fo0, as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
  mfExt <- data.frame(1,U,Vxb)
  colnames(mfExt)<-c("(Intercept)",nnomi)
  objF <- lm(Fo, data = mfExt)
  
  #browser()
  
  objW<-objF
  
  #controllo se qualche coeff e' NA..
  isNAcoef<-any(is.na(objF$coefficients))
 
  #browser()
  if(isNAcoef) {
    nameNA.psi <- names(objF$coefficients)[which(is.na(objF$coefficients))]
    nameNA.U <- gsub("psi", "U", nameNA.psi)
    if(fix.npsi) {
        cat("breakpoint estimate(s):", as.vector(psi), "\n")
        stop("coef ", nameNA.psi, " is NA: breakpoint(s) at the boundary or too close together", call. = FALSE)
    } else {
        warning("some estimate is NA (too many breakpoints?): removing ",
            length(nameNA.psi), " jump-point(s)", call. = FALSE)
        Fo <- update(Fo, paste(".~ .-", nameNA.U, "-", nameNA.psi))
        objF <- lm(Fo, data = mfExt) 
        idNA.psi <- match(nameNA.psi, nomiVxb)
        nomiVxb <- setdiff(nomiVxb, nameNA.psi)
        nomiU <- setdiff(nomiU, nameNA.U)
        Xtrue <- Xtrue[, -idNA.psi, drop = FALSE]
        PSI.mid<- PSI.mid[, -idNA.psi, drop = FALSE]
        id.psi.group <- id.psi.group[-idNA.psi]
        psi <- psi[-idNA.psi]
        psi.rounded <- psi.rounded[, -idNA.psi, drop = FALSE]
    }
  }

  #organizziamo i risultati da restituire per psi..
  colnames(psi.rounded)<-names(psi)<-nomiVxb
  rownames(psi.rounded)<-c("inf [","sup (")
  # Cov <- vcov(objF) 
  # 
  # var.Tay<-function(est1,est2,v1,v2,v12){
  #   r<- est1/est2
  #   vv<-(v1+v2*r^2-2*r*v12)/est2^2
  #   vv}
  # 
  # 
  # #browser()
  # 
  # #var.Tay(num, den, v.g, v.b, cov.g.b)
  # varPsi<- rep(NA, length(nomiU))
  # for(j in 1:length(nomiU)){
  #   num<-objF$coefficients[nomiVxb[j]]
  #   den<-objF$coefficients[nomiU[j]]
  #   v.g <-Cov[nomiVxb[j],nomiVxb[j]]
  #   v.b<- Cov[nomiU[j],nomiU[j]]
  #   cov.g.b <- Cov[nomiVxb[j],nomiU[j]]
  #   #if(is.null(rho)) {
  #     rho<-mean(Xtrue[, nomiZ[j] ,drop=TRUE]<psi[[nomiVxb[j]]])
  #     #rho<- 1-exp(-rho*(n^(1/3))) #rho^(2*sqrt(2/n)) #1-exp(-5*rho) #con 5,6 valori piu' alti => SE piu' piccoli..
  #     rho<-  rho^(sqrt(1/n))
  #   #}
  #   cov.g.b<- rho*sqrt(v.g*v.b)
  #   varPsi[j]<-var.Tay(num, den, v.g, v.b, cov.g.b)
  # }
  # names(varPsi) <- nomiVxb
  # 
  # #browser()
  # Cov[nomiVxb, ]<- Cov[, nomiVxb] <- 0
  # diag(Cov)[nomiVxb]<-varPsi
  # #Cov[nomiVxb, nomiVxb ]<- varPsi
  # 
  # 
  # #browser()
  # #var.Tay(num, den, v.g, v.b, cov.g.b)
  # 
  # id <- match(nomiVxb, names(coef(objF)))
  # vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
  ris.psi <-matrix(NA,length(psi),3)
  colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
  rownames(ris.psi) <- nomiVxb
  ris.psi[,2]<-psi
  #ris.psi[,3]<-sqrt(vv)
  
  ##  solo per simulazioni
  #browser()
  
#  ris.psi<-cbind(ris.psi,
 #   st0=sqrt(var.Tay(num, den, v.g, v.b, 0)),
  #  st99=sqrt(var.Tay(num, den, v.g, v.b, .99*sqrt(v.g*v.b))))
  
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
  #browser()
  id.psi<- x%in%psi.rounded[1,]
  #=================================================
  ##RI-AGGIUNGI IL MINIMO!!!!!!!!!!
  if(adjX){ #ATTENZIONE.. e se ci sono piu' breakpoints o piu' variabili (con piu' breakpoints)??
    psi.rounded<- psi.rounded + min.x
    ris.psi[,2] <- ris.psi[,2] + min.x
  }
  
  objF$psi <- ris.psi

  a<-rep(1:Tsp[3L], l=length(y))
  b<-rep(Tsp[1L]:Tsp[2L], each=Tsp[3L])[seq_len(length(y))]
  break.dates <- paste(b,"(",a,")",sep="")[id.psi]
  attr(psi.rounded,"break.dates") <- break.dates
  objF$psi.rounded <- psi.rounded
  
  #stima il modello "vero" (non-working)
  U <- (Xtrue > PSI.mid)
  colnames(U)<-nomiU
  X <- cbind(x.lin,U)
  objF$obj.ok <- mylm(X, y) #coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b))
  objF$objW<- objW
  objF$fitted.values<-objF$obj.ok$fitted.values
  objF$residuals<- objF$obj.ok$residuals
  objF$coefficients[1:length(objF$obj.ok$coefficients)] <- objF$obj.ok$coefficients
  objF$coefficients[nomiVxb] <-psi.rounded[1,]
  objF$nameUV <- list(U = drop(nomiU), V = nomiV, Z = name.Z) #Z = name.Z
  objF$rangeZ<-obj$rangeZ
  objF$Z <- Z[,unique(name.Z),drop=FALSE]
  if(n.boot>0) objF$seed <- seed
  if(adjX) {
    objF$Z <- objF$Z + min.x
    objF$rangeZ<- objF$rangeZ + min.x
  }

  objF$call <- match.call()
  objF$orig.call<-orig.call
  objF$psi.history <- psi.values
  objF$it <- it 
  objF$epsilon <- obj$epsilon
  objF$id.warn <- id.warn
  #objF$rho<-rho
  objF$psi<- objF$psi[,-1,drop=FALSE] #rimuovi la colonna Initial
  if(var.psi){
    Cov <- vcov.stepmented(objF, k=NULL)
    id <- match(nomiVxb, names(coef(objF)))
    vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
    objF$psi[,"St.Err"]<-sqrt(vv)
    objF$vcov<- Cov
  }
  class(objF) <- c("stepmented","lm")
  return(objF)
}
