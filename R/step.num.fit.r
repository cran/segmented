step.num.fit<-function(y, x.lin, Xtrue, PSI, ww, opz, return.all.sol=FALSE){  
  #----------------------
  search.min<-function(h, psi, psi.old, X, y, w, id.fix.psi=NULL) {
    psi.ok<- psi*h + psi.old*(1-h)
    psi.ok[id.fix.psi]<- psi.old[id.fix.psi]
    PSI <- matrix(rep(psi.ok, rep(n, length(psi.ok))), ncol = length(psi.ok))
    U1 <- (Xtrue>PSI) #(Z - PSI) * (Z > PSI)
    #if (pow[1] != 1) U1 <- U1^pow[1]
    obj1 <- try(mylm(cbind(X, U1), y, w), silent = TRUE)
    if (class(obj1)[1] == "try-error") obj1 <- try(lm.wfit(cbind(X, U1), y, w), silent = TRUE)
    #if (class(obj1)[1] == "try-error") obj1 <- try(.lm.fit(cbind(X, U1), y), silent = TRUE)
    L1 <- if (class(obj1)[1] == "try-error") L0 + 10
    else sum(w*obj1$residuals^2)
    #r<-sum(obj1$residuals^2 * w)
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
  mylm <-function(x,y,w=1){ #,w=1,offs=0 in step.st.fit() non ci sovrebbero essere w e offs
    x1<-x*sqrt(w)
    #y<-y-offs
    y1<-y*sqrt(w)
    #o<-.lm.fit(y=y,x=x)
    #b<-o$coefficients
    #fit<- o$fitted.values
    #r<-o$residuals
    b<-drop(solve(crossprod(x1),crossprod(x1,y1))) #.lm.fit(x=x1, y=y1) #lm.wfit(x, y, w) #CHI E' PIU' VELOCE?
    #b<-solve(crossprod(x),crossprod(x,y))
    #browser()
    #fit<- drop(tcrossprod(x,t(b)))
    fit<- drop(tcrossprod(b,x))
    r<-y-fit
    o<-list(coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b))
    o
  }
  #-----------
  adj.psi <- function(psii, LIM) {
    pmin(pmax(LIM[1, ], psii), LIM[2, ])
  }
  
  #------------
  tol<-opz$toll
  display<-opz$display
  it.max<-opz$it.max
  dev0<-opz$dev0
  useExp.k<-opz$useExp.k
  min.step<- opz$min.step #=.0001
  conv.psi<-opz$conv.psi #=FALSE
  alpha<-opz$alpha
  limZ <- apply(Xtrue, 2, quantile, names = FALSE, probs = c(alpha[1], alpha[2]))
  #limZ <- apply(Xtrue, 2, quantile, names = FALSE, probs = c(alpha, 1 - alpha))
  fix.npsi<-opz$fix.npsi
  agg<-opz$agg
  h<-opz$h
  npsii<-opz$npsii
  npsi<- sum(npsii) #opz$npsi
  P<-length(npsii) #P<-opz$P
  digits <- opz$digits
  rangeZ <- opz$rangeZ
  
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
  
  
  if(it.max==0){
    obj <- lm.wfit(x = cbind(x.lin, Xtrue>PSI), y = y, w = ww)
    L1 <- sum(obj$residuals^2 * ww)
    obj$epsilon <- epsilon
    idZ<-(plin+1):(plin+ncol(PSI))
    b<- obj$coef[idZ]
    obj <- list(obj = obj, psi = PSI[1,], psi.values = psi.values, 
                rangeZ = rangeZ, beta.c=b, epsilon = epsilon,  
                SumSquares.no.gap = L1,  
                id.warn = TRUE)
    return(obj)
  }
  
  #PSI0<- matrix(psi0, n, npsi, byrow = TRUE)
  #XREG <- cbind(x.lin, Xtrue>PSI)
  #obj0 <- try(mylm(XREG, y), silent = TRUE) 
  #L0<- sum(obj0$residuals^2) #*ww
  L0<- dev0*.8 #e' inutil efar stimare un modello per avere L0 (che eserve solo se display=TRUE). Lo stesso in step.ts.fit
  n.intDev0<-nchar(strsplit(as.character(L0),"\\.")[[1]][1])
  dev.values[length(dev.values) + 1] <- dev0#opz$dev0 #del modello iniziale (senza psi)
  dev.values[length(dev.values) + 1] <- L0 #modello con psi iniziali
  psi0<-PSI[1,]
  psi.values[[length(psi.values) + 1]] <- psi0 #psi iniziali
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
  low <- apply(Xtrue, 2, min)
  up <- apply(Xtrue, 2, max)

  L1<-L0+10  
  #==============================================
  while (abs(epsilon) > tol) {
    i <- i + 1
    #if(i==1) browser()
    xx <- Xtrue[,cumsum(npsii),drop=FALSE]
    for (p in 1:P) {
      psis <- sort(psi0[pos[[p]]])
      gruppi <- cut(xx[,p], breaks = c(low[p] - 0.1, psis, up[p]), labels = FALSE)
      points <- c(low[p], psis, up[p])
      right <- c(low[p], points[2:(npsii[p] + 1)] + agg[pos[[p]]][order(psi0[pos[[p]]])] * (points[3:(npsii[p] + 2)] - points[2:(npsii[p] + 1)]), NA)
      left <- c(NA, points[2:(npsii[p] + 1)] - agg[pos[[p]]][order(psi0[pos[[p]]])] * (points[2:(npsii[p] + 1)] - points[1:npsii[p]]), up[p])
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
    
    obj <- lm.wfit(y = y, x = XREG, w=ww)
    
    #b <- obj$coef[(2:(sum(k) + 1))]
    #g <- obj$coef[((sum(k) + 2):(2 * sum(k) + 1))]
    
    idZ<-(plin+1):(plin+ncol(Z))
    idW<-(plin+ncol(Z)+1): ( plin+ncol(Z)+ncol(W))
    b<- obj$coef[idZ]
    g<- obj$coef[idW]
    
    if(any(is.na(c(b, g)))){
      if(return.all.sol) return(list(dev.values, psi.values)) else stop("breakpoint estimate too close or at the boundary causing NA estimates.. too many breakpoints being estimated?", call.=FALSE)
    }
    psi1 <- -g/b
    psi1<- psi0+ opz$h*(psi1-psi0)
    psi1<- adj.psi(psi1, limZ) #limZ rangeZ??
    psi1<-unlist(tapply(psi1, opz$id.psi.group, sort), use.names =FALSE)
    
    #if(i==2) browser()
    #la f e' chiaramente a gradino per cui meglio dividere..
    a0<-optimize(search.min, c(0,.5), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww)
    a1<-optimize(search.min, c(.5,1), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww)
    a<-if(a0$objective<=a1$objective) a0 else a1
    
    #M<-1
    #while(L1>L0){
    #  a<-optimize(search.min, c(0,M), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, offs=offs)
    #  L1<- a$objective
    #  M<-M*.3
    #}
    if(a$objective<L0){
      k.values[length(k.values) + 1] <- use.k <- a$minimum
      L1<- a$objective
    } else {
      k.values[length(k.values) + 1] <- use.k <- 0
      L1<- L0  
    }
    if(use.k<=.01){
      k.List<-j.List<-NULL
      for(j in 1:length(psi1)){
        a0<-optimize(search.min, c(0,.5), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, id.fix.psi=j)
        a1<-optimize(search.min, c(.5,1), psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww, id.fix.psi=j)
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
        L1=search.min(1, psi=psi1, psi.old=psi0, X=x.lin, y=y, w=ww)
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
                "  est.psi = ",paste(formatC(unlpsi,digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                sep=""), "\n")
    }
    
    epsilon <- if(conv.psi) max(abs((psi1 -psi0)/psi0)) else (L0 - L1)/(abs(L0) + 0.1) 
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
  
  #ATTENZIONE .. Assume che obj sia stato stimato sempre!
  obj<-list(obj=obj, psi=psi1, psi.values=psi.values, rangeZ=rangeZ, SumSquares.no.gap=L1, 
            beta.c=b, it=i, epsilon=epsilon, id.warn=id.warn) 
  return(obj)
} #end jump.fit
