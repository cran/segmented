step.glm.fit<-function(y, x.lin, Xtrue, PSI, ww, offs, opz, return.all.sol=FALSE){  
  #----------------------
  search.min<-function(h, psi, psi.old, X, y, w, offs, id.fix.psi=NULL) {
    psi.ok<- psi*h + psi.old*(1-h)
    psi.ok[id.fix.psi]<- psi.old[id.fix.psi]
    PSI <- matrix(rep(psi.ok, rep(n, length(psi.ok))), ncol = length(psi.ok))
    U1 <- (Xtrue>PSI) #(Z - PSI) * (Z > PSI)
    #if (pow[1] != 1) U1 <- U1^pow[1]
    obj1 <- try(suppressWarnings(glm.fit(x = cbind(X, U1), y = y, offset = offs,
                      weights = w, family = fam, control = glm.control(maxit = maxit.glm1[i]), etastart = eta0)),
                      silent = TRUE)
    L1 <- if (class(obj1)[1] == "try-error") L0 + 10 else obj1$dev
    attr(L1, "eta") <- obj1$linear.predictor
    L1
  }
  toMatrix<-function(.x, ki){
    # ripete ogni .x[,j] ki[j] volte
    if(ncol(.x)!=length(ki)) stop("It should be ncol(.x)==length(ki)")
    if(all(ki==1)) return(.x)
    M<-vector("list", length=length(ki))
    for(j in 1:length(ki)) M[[j]]<-replicate(ki[[j]], cbind(.x[,j]), simplify=TRUE)
    do.call(cbind, M)
  }
  ### -----
  # mylm<-function(x,y,w=1,offs=0){
  #   x1<-x*sqrt(w)
  #   y<-y-offs
  #   y1<-y*sqrt(w)
  #   b<-drop(solve(crossprod(x1),crossprod(x1,y1)))
  #   fit<-drop(tcrossprod(x,t(b)))
  #   r<-y-fit
  #   o<-list(coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b))
  #   o
  # }
  #-----------
  adj.psi <- function(psii, LIM) {
    pmin(pmax(LIM[1, ], psii), LIM[2, ])
  }
  #------------
  #-----------
  
  fam<-opz$fam
  maxit.glm<-opz$maxit.glm
  #--------------
  tol<-opz$toll
  display<-opz$display
  it.max<-opz$it.max
  #dev0<-opz$dev0
  useExp.k<-opz$useExp.k
  min.step<- opz$min.step #=.0001
  conv.psi<-opz$conv.psi #=FALSE
  alpha<-opz$alpha
  #limZ <- apply(Xtrue, 2, quantile, names = FALSE, probs = c(alpha[1], alpha[2]))
  limZ <- if(is.null(opz$limZ)) apply(Xtrue, 2, quantile, names=FALSE, probs=alpha) else opz$limZ
  
  fix.npsi<-opz$fix.npsi
  agg<-opz$agg
  hh <-opz$h
  npsii<-opz$npsii
  npsi<- sum(npsii) #opz$npsi
  P<-length(npsii) #P<-opz$P
  digits<-opz$digits
  rangeZ<-opz$rangeZ
  
  #  pos.vec <- 1:npsi
  #  pos <- vector("list", P)
  #  ind <- 0
  pos<- tapply(1:npsi, rep(1:P, npsii), list)
  i <- 0
  agg <- rep(agg, npsi)
  #  direz <- matrix(NA, it.max, npsi)
  #  conv <- rep(FALSE, npsi)
  #  ind.conv <- NULL
  n<-length(y)
  plin<-ncol(x.lin)
  epsilon<-10
  k.values<-dev.values<- NULL
  psi.values <-list()
  psi.values[[length(psi.values) + 1]] <- NA
  #PSI0<- matrix(psi0, n, npsi, byrow = TRUE)
  XREG <- cbind(x.lin, Xtrue>PSI)

  if(it.max==0){
    obj <- suppressWarnings(glm.fit(x = XREG, y = y, offset = offs,
                                    weights = ww, family = fam))
    L1 <- obj$dev 
    obj$epsilon <- epsilon
    idZ<-(plin+1):(plin+ncol(PSI))
    b<- obj$coef[idZ]
    obj <- list(obj = obj, psi = PSI[1,], psi.values = psi.values, 
                rangeZ = rangeZ, beta.c=b, epsilon = epsilon,  
                SumSquares.no.gap = L1,  
                id.warn = TRUE)
    return(obj)
  } 
  
  if(!opz$usestepreg){
    dev.values[length(dev.values) + 1] <- opz$dev0 #modello senza psi 
    psi.values[[length(psi.values) + 1]] <- NA #nessun psi 
  }
  
  if(is.null(opz$fit.psi0)){
    obj <- suppressWarnings(glm.fit(x = XREG, y = y, offset = offs,
                                    weights = ww, family = fam, etastart=opz$eta0))
    L0 <- obj$dev
    eta0 <- obj$linear.predictors
  } else {
    L0   <- opz$fit.psi0$L0
    eta0 <- opz$fit.psi0$eta0
  }
    

  n.intDev0<-nchar(strsplit(as.character(L0),"\\.")[[1]][1])
  #dev.values[length(dev.values) + 1] <- dev0#opz$dev0 #del modello iniziale (senza psi)
  dev.values[length(dev.values) + 1] <- L0 #modello con psi iniziali
  psi0<-PSI[1,]
  psi.values[[length(psi.values) + 1]] <- psi0 #psi iniziali
  
  if(is.null(maxit.glm)){
    Nboot <- if(is.null(opz$Nboot)) 0 else opz$Nboot
    maxit.glm1 <- rep(1:it.max + Nboot, 1:it.max+1) #2*rep(1:it.max, 1:it.max)
    maxit.glm1 <- pmin(maxit.glm1, 25)
  } else {
    maxit.glm1 <- rep(maxit.glm, it.max)
  }
  
  
  #==============================================
  if (display) {
    unlpsi<- unlist(psi0)
    Lp<-length(unlpsi)
    
    cat(paste("iter = ", sprintf("%2.0f",0),
              "  dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), L0), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
              "  k = ", sprintf("%5.0f", NA),
              "  n.psi = ",formatC(Lp,digits=0,format="f"), 
              "  ini.psi = ",paste(formatC(unlpsi[1:min(5,Lp)],digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
              sep=""), "\n")
  }
  id.warn <- FALSE
  
  low <- apply(Xtrue[,unique(colnames(Xtrue)),drop=FALSE], 2, min)
  up <-  apply(Xtrue[,unique(colnames(Xtrue)),drop=FALSE], 2, max)
  
  
  L1<-L0+10
  tolOp<-if(is.null(opz$tol.opt)) seq(.001, .Machine$double.eps^0.25, l=it.max) else rep(opz$tol.opt, it.max)
  #==============================================
  
  idZ<-(plin+1):(plin+ncol(PSI))
  idW<-(plin+ncol(PSI)+1): ( plin+2*ncol(PSI))
  
  while (abs(epsilon) > tol) {
    i <- i + 1
    #if(i==1) browser()
    xx <- Xtrue[,cumsum(npsii),drop=FALSE]
    for (p in 1:P) {
      psis <- sort(psi0[pos[[p]]])
      gruppi <- cut(xx[,p], breaks = c(low[p] - 0.1, psis, up[p]), labels = FALSE)
      if(any(is.na(gruppi))) stop(paste("too many breaks for step term #", p, "?"), call.=TRUE)
      points <- c(low[p], psis, up[p])
      right <- c(low[p], points[2:(npsii[p] + 1)] + agg[pos[[p]]][order(psi0[pos[[p]]])] * (points[3:(npsii[p] + 2)] - points[2:(npsii[p] + 1)]), NA)
      left <- c(NA, points[2:(npsii[p] + 1)] - agg[pos[[p]]][order(psi0[pos[[p]]])] * (points[2:(npsii[p] + 1)] - points[1:npsii[p]]), up[p])
      #if(any(is.na(left))| any(is.na(right))) stop(paste("too many breaks for step term #", p, "?"), call.=TRUE)
      for (j in 1:(npsii[p] + 1)) {
        xx.j <- xx[,p][gruppi == j]
        xx[,p][gruppi == j] <- right[j] + (xx.j - points[j]) * 
          ((left[j + 1] - right[j])/(points[j + 1] - points[j]))
      }
    }
    
    XX<-toMatrix(xx, npsii)
    PSI<- matrix(psi0, n, npsi, byrow = TRUE)
    W <- (1/(2 * abs(XX - PSI)))
    Z <- (XX * W + 1/2)
    XREG <- cbind(x.lin, Z, W)
    
    #obj<-try(mylm(XREG,y,w=ww,offs=offs), silent = TRUE)
    #if(class(obj)[1]=="try-error") 
    #  obj <- lm.wfit(y = y, x = XREG, offset = offs, w=ww )
    #b <- obj$coef[(2:(sum(k) + 1))]
    #g <- obj$coef[((sum(k) + 2):(2 * sum(k) + 1))]
    obj <- suppressWarnings(glm.fit(x = XREG, y = y, offset = offs,
                                     weights = ww, family = fam, control = glm.control(maxit = maxit.glm1[i]), etastart = eta0)) 
      
      
    
    #idZ<-(plin+1):(plin+ncol(Z))
    #idW<-(plin+ncol(Z)+1): ( plin+ncol(Z)+ncol(W))
    b<- obj$coef[idZ]
    g<- obj$coef[idW]
    
    if(any(is.na(c(b, g)))){
      if(return.all.sol) return(list(dev.values, psi.values)) else stop("breakpoint estimate too close or at the boundary causing NA estimates.. too many breakpoints being estimated?", call.=FALSE)
    }
    
    psi1 <- -g/b
    psi1<- psi0+ opz$h*(psi1-psi0)
    #aggiusta la stima di psi..
    psi1<- adj.psi(psi1, limZ)
    psi1<-unlist(tapply(psi1, opz$id.psi.group, sort), use.names =FALSE)
    
    #if(i==1) browser()
    #la f e' chiaramente a gradino per cui meglio dividere..
     a0<-optimize(search.min, c(0,.5), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs, tol=tolOp[i])
     a1<-optimize(search.min, c(.5,1), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs, tol=tolOp[i])
     a <-if(a0$objective<=a1$objective) a0 else a1
     
     
     
     #a0<-optimize(search.min, c(0,.33), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs)
     #a1<-optimize(search.min, c(.33,.66), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs)
     #a2<-optimize(search.min, c(.66,1), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs)
     #a<-if(a2$objective<=a$objective) a2 else a
    
    if(a$objective<L0){
      k.values[length(k.values) + 1] <- use.k <- a$minimum
      L1<- a$objective
      eta0<- attr(a$objective, "eta")
    } else {
      k.values[length(k.values) + 1] <- use.k <- 0
      L1<- L0  
    }
     
    if(use.k<=.01){
      k.List<-j.List<-NULL
     for(j in 1:length(psi1)){
       a0<-optimize(search.min, c(0,.5), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs, id.fix.psi=j, tol=tolOp[i])
       a1<-optimize(search.min, c(.5,1), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs, id.fix.psi=j, tol=tolOp[i])
       a <-if(a0$objective<=a1$objective) a0 else a1
       if(a$objective<L1){
         j.List[[j]]<-setdiff(1:length(psi1),j) #indici di psi che devono cambiare..
         k.List[[j]]<-a$minimum
       } else {
         j.List[[j]]<-NA
         k.List[[j]]<-NA
         
       }
     }
      id.to.be.changed<- unique(unlist(j.List[!sapply(k.List, is.na)]))
      if(!is.null(id.to.be.changed)){
        use.k<-rep(0,length(psi1))
        use.k[id.to.be.changed] <-mean(unlist(k.List[!sapply(k.List, is.na)]))
        psi1 <- psi1*use.k + psi0* (1-use.k)
        use.k<-mean(use.k)
        L1=search.min(1, psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs)
      } else {
        psi1<-psi0
      }
    } else {
      psi1 <- psi1*use.k + psi0* (1-use.k)
    }
     

    if (!is.null(digits)) psi1 <- round(psi1, digits)
    #PSI1 <- matrix(psi1, n, npsi, byrow = TRUE)
    #XREG1 <- cbind(x.lin, Xtrue>PSI1)
    #obj1 <- try(mylm(XREG1, y, ww, offs), silent = TRUE)
    #if (class(obj1)[1] == "try-error") obj1 <- try(lm.wfit(XREG1, y, ww, offs), silent = TRUE)
    delta<- psi1-psi0
    
    if (display) {
      flush.console()
      #n.intDev0<-nchar(strsplit(as.character(dev.values[2]),"\\.")[[1]][1])
      unlpsi<- unlist(psi1)
      Lp<-length(unlpsi)
      
      cat(paste("iter = ", sprintf("%2.0f",i),
                "  dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), L1), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
                "  k = ", sprintf("%2.3f", use.k),
                "  n.psi = ",formatC(Lp,digits=0,format="f"), 
                "  est.psi = ",paste(formatC(unlpsi[1:min(Lp,5)],digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                sep=""), "\n")
    }
    
    epsilon <- (L0 - L1)/(abs(L0) + 0.1) 
    L0<-L1
    
    k.values[length(k.values)+1]<-use.k
    psi.values[[length(psi.values) + 1]] <- psi1
    dev.values[length(dev.values) + 1] <- L0
    
    if (i >= it.max) {
      id.warn <- TRUE
      break
    }
    psi0<-psi1
  } #end while_it
  
  #browser()
  psi1 <-unlist(tapply(psi1, opz$id.psi.group, sort))
  PSI<- matrix(psi1, n, npsi, byrow = TRUE)
  U <- 1*(Xtrue>PSI)
  
  #ATTENZIONE .. Assume che obj sia stato stimato sempre!
  obj<-list(obj=obj, psi=psi1, psi.values=psi.values, rangeZ=rangeZ, SumSquares.no.gap=L1, 
            beta.c=b, it=i, epsilon=epsilon, id.warn=id.warn, U=U, eta0=eta0) 
  return(obj)
} #end jump.fit
