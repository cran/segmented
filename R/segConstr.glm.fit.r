segConstr.glm.fit <-function (y, XREG, Z, PSI, w, offs, opz, return.all.sol = FALSE) 
{
    useExp.k = TRUE
    search.min <- function(h, psi, psi.old, X, y, w, offs) {#DUBBIO: Ma fam, eta0, L0 e maxit.glm devo passarli come argom,enti o li trova?
      psi.ok<- psi*h + psi.old*(1-h)
      PSI <- matrix(rep(psi.ok, rep(n, length(psi.ok))), ncol = length(psi.ok))
      U1 <- (Z - PSI) * (Z > PSI)
      #if (pow[1] != 1) U1 <- U1^pow[1]
      for(i in 1:length(RList)){#trasforma le U
        UList[[i]]<- cbind(Zseg[,i], U1[, id.psi.group==i])%*%invA.RList[[i]] #
        #nomiUList[[i]]<- rep(i, ncol(UList[[i]]) )
      }
      U1<-do.call(cbind, UList) #la matrice del disegno sara' cbind(X, U1)
      
      obj1 <- try(suppressWarnings(glm.fit(x = cbind(X, U1), y = y, offset = offs,
                weights = w, family = fam, control = glm.control(maxit = maxit.glm), etastart = eta0)),
                silent = TRUE)
      L1 <- if (class(obj1)[1] == "try-error") L0 + 10 else obj1$dev
      L1
    }
    est.k<-function(x1,y1,L0){
      ax<-log(x1)
      .x<-cbind(1,ax,ax^2)
      b<-drop(solve(crossprod(.x),crossprod(.x,y1)))
      const<-b[1]-L0
      DD<-sqrt(b[2]^2-4*const*b[3])
      kk<-exp((-b[2]+ DD) /(2*b[3]))
      return(round(kk))
      
      #  ff<-function(xx) b[1]+b[2]*xx + b[3]*xx^2+ L0
      #  a<-uniroot(ff, c(log(x[4]), 3.4))
    }
    dpmax <- function(x, y, pow = 1) {
        if (pow == 1) 
            -(x > y)
        else -pow * ((x - y) * (x > y))^(pow - 1)
    }
    in.psi <- function(LIM, PSI, ret.id = TRUE) {
        a <- PSI[1, ] <= LIM[1, ]
        b <- PSI[1, ] >= LIM[2, ]
        is.ok <- !a & !b
        if (ret.id) 
            return(is.ok)
        isOK <- all(is.ok) && all(!is.na(is.ok))
        isOK
    }
    far.psi <- function(Z, PSI, id.psi.group, ret.id = TRUE, fc = 0.93) {
        nSeg <- length(unique(id.psi.group))
        npsij <- tapply(id.psi.group, id.psi.group, length)
        nj <- sapply(unique(id.psi.group), function(.x) {
            tabulate(rowSums((Z > PSI)[, id.psi.group == .x, 
                drop = FALSE]) + 1)
        }, simplify = FALSE)
        ff <- id.far.ok <- vector("list", length = nSeg)
        for (i in 1:nSeg) {
            if (length(nj[[i]]) != npsij[i] + 1) 
                nj[[i]] <- tabulate(rowSums((Z >= PSI)[, id.psi.group == 
                  i, drop = FALSE]) + 1)
            id.ok <- (nj[[i]] >= 2)
            id.far.ok[[i]] <- id.ok[-length(id.ok)] & id.ok[-1]
            ff[[i]] <- ifelse(diff(nj[[i]]) > 0, 1/fc, fc)
        }
        id.far.ok <- unlist(id.far.ok)
        ff <- unlist(ff)
        if (!ret.id) {
            return(all(id.far.ok))
        }
        else {
            attr(id.far.ok, "factor") <- ff
            return(id.far.ok)
        }
    }
    adj.psi <- function(psii, LIM) {
        pmin(pmax(LIM[1, ], psii), LIM[2, ])
    }
    #nuovo per i vincoli

    RList <- opz$RList
    nomiUList<-UList<- vector("list",length(RList))
    invAList <- lapply(RList, function(.x)rbind(c(1,rep(0,nrow(.x)-1)),diff(diag(nrow(.x)))))
    invA.RList<-lapply(1:length(RList), function(i) invAList[[i]]%*% RList[[i]])
    nomiUList<- lapply(1:length(RList), function(i)rep(i, ncol(RList[[i]])))

    #-----------
    eta0<-opz$eta0
    fam<-opz$fam
    maxit.glm<-opz$maxit.glm
    #----------------------------
    
    n<-length(y)
    min.step<-opz$min.step
    rangeZ <- apply(Z, 2, range)
    alpha<-opz$alpha
    limZ <- apply(Z, 2, quantile, names=FALSE, probs=c(alpha[1],alpha[2]))
    psi<-PSI[1,]
    id.psi.group<-opz$id.psi.group
    conv.psi<-opz$conv.psi 
    digits<-opz$digits
    pow<-opz$pow
    nomiOK<-opz$nomiOK
    toll<-opz$toll
    hh<-opz$h
    gap<-opz$gap
    #fix.npsi<-opz$fix.npsi
    fix.npsi<-opz$stop.if.error
    dev.new<-opz$dev0
    visual<-opz$visual
    it.max<-old.it.max<-opz$it.max
    fc<-opz$fc
    names(psi)<-id.psi.group
    it <- 0
    epsilon <- 10
    k.values<-dev.values<- NULL
    psi.values <-list()
    psi.values[[length(psi.values) + 1]] <- NA
    #id.psi.ok<-rep(TRUE, length(psi))
    sel.col.XREG<-unique(sapply(colnames(XREG), function(x)match(x,colnames(XREG))))
    if(is.numeric(sel.col.XREG)) XREG<-XREG[,sel.col.XREG,drop=FALSE] #elimina le ripetizioni, ad es. le due intercette..
    #invXtX <- opz$invXtX
    #Xty <- opz$Xty
    if (!in.psi(limZ, PSI, FALSE)) 
        stop("starting psi out of the range.. see 'alpha' in seg.control.", 
            call. = FALSE)
    if (!far.psi(Z, PSI, id.psi.group, FALSE)) 
        stop("psi values too close each other. Please change (decreases number of) starting values", 
            call. = FALSE)
    n.psi1 <- ncol(Z)
    
    #browser()
    
    Zseg <- XREG[,opz$nomiSeg,drop=FALSE] #
    XREG <- XREG[, -match(opz$nomiSeg, colnames(XREG)),drop=FALSE]
    
    U <- ((Z - PSI) * (Z > PSI))
    #if (pow[1] != 1) U <- U^pow[1]
    
    for(i in 1:length(RList)){#trasforma le U
      UList[[i]]<- cbind(Zseg[,i], U[, id.psi.group==i])%*%invA.RList[[i]] #
      #nomiUList[[i]]<- rep(i, ncol(UList[[i]]) )
    }
    U<-do.call(cbind, UList) #la matrice del disegno sara' cbind(X, U)
    
    obj0 <- suppressWarnings(glm.fit(x = cbind(XREG, U), y = y, offset = offs,
                                     weights = w, family = fam, control = glm.control(maxit = maxit.glm), etastart = eta0))
    eta0<- obj0$linear.predictors
    L0<- obj0$dev
    n.intDev0<-nchar(strsplit(as.character(L0),"\\.")[[1]][1])
    dev.values[length(dev.values) + 1] <- opz$dev0 #del modello iniziale (senza psi)
    dev.values[length(dev.values) + 1] <- L0 #modello con psi iniziali
    psi.values[[length(psi.values) + 1]] <- psi #psi iniziali

    if (visual) { #questo e' il visual di "lm"
        cat(paste("iter = ", sprintf("%2.0f", 0), "  dev = ", 
            sprintf(paste("%", n.intDev0 + 6, ".5f", sep = ""), 
                L0), "  k = ", sprintf("%2.0f", NA), "  n.psi = ", 
            formatC(length(unlist(psi)), digits = 0, format = "f"), 
            "  ini.psi = ", paste(formatC(unlist(psi), digits = 3, 
                format = "f"), collapse = "  "), sep = ""), "\n")
    }
    id.warn <- FALSE
    id.psi.changed <- rep(FALSE, it.max)
    #============================================== inizio ciclo
    #browser()
    #Zseg <- XREG[,opz$nomiSeg,drop=FALSE]
    #XREG <- XREG[, -match(opz$nomiSeg, colnames(XREG)),drop=FALSE]
    while (abs(epsilon) > toll) {
      it<-it+1
      n.psi0 <- n.psi1
      n.psi1 <- ncol(Z)
      if(n.psi1!=n.psi0){
        U <- ((Z-PSI)*(Z>PSI)) #pmax((Z - PSI), 0)^pow[1]
        if(pow[1]!=1) U<-U^pow[1]
        obj0 <- suppressWarnings(glm.fit(x = cbind(XREG, U), y = y, offset = offs,
                                         weights = w, family = fam, control = glm.control(maxit = maxit.glm), etastart = eta0))
        eta0<-obj0$linear.predictors
        L0< - obj0$dev
      } else {
        V <- (Z>PSI)
        U <- (Z - PSI) * V
        V <- -V
      }
      #V <- dpmax(Z,PSI,pow=pow[2])# ifelse((Z > PSI), -1, 0)
      
        for(i in 1:length(RList)){#trasforma le U
          UList[[i]]<- cbind(Zseg[,i], U[, id.psi.group==i])%*%invA.RList[[i]] #
          nomiUList[[i]]<- rep(i, ncol(UList[[i]]) )
        }
        U<-do.call(cbind, UList)

        X <- cbind(XREG, U, V)
        rownames(X) <- NULL
        #colnames(X)[(ncol(XREG) + 1):ncol(U)] <- paste("U", 
         #   1:ncol(U), sep = "") #, paste("V", 1:ncol(V), sep = ""))
        obj <- suppressWarnings(glm.fit(X, y, offset = offs, weights = w, family = fam, 
                                        control = glm.control(maxit = maxit.glm), etastart = eta0))
        eta0<-obj$linear.predictors
        beta.c <- coef(obj)[ncol(XREG)+(1:ncol(U))]
        coefUList <- lapply(1:length(RList), function(i) (invA.RList[[i]]%*%beta.c[unlist(nomiUList)==i])[-1])
        beta.c <- unlist(coefUList)

        gamma.c <- coef(obj)[colnames(Z)]
        if (any(is.na(c(beta.c, gamma.c)))) {
            if (fix.npsi) {
                if (return.all.sol) 
                  return(list(dev.values, psi.values))
                else stop("breakpoint estimate too close or at the boundary causing NA estimates.. too many breakpoints being estimated?", 
                  call. = FALSE)
            } else {
                id.coef.ok <- !is.na(gamma.c)
                psi <- psi[id.coef.ok]
                if (length(psi) <= 0) {
                  warning(paste("All breakpoints have been removed after", 
                    it, "iterations.. returning 0"), call. = FALSE)
                  return(0)
                }
                gamma.c <- gamma.c[id.coef.ok]
                beta.c <- beta.c[id.coef.ok]
                Z <- Z[, id.coef.ok, drop = FALSE]
                rangeZ <- rangeZ[, id.coef.ok, drop = FALSE]
                limZ <- limZ[, id.coef.ok, drop = FALSE]
                nomiOK <- nomiOK[id.coef.ok]
                id.psi.group <- id.psi.group[id.coef.ok]
                names(psi) <- id.psi.group
            }
        }
        psi.old <- psi
        psi <- psi.old + hh*gamma.c/beta.c
        #aggiusta la stima di psi..
        psi<- adj.psi(psi, rangeZ)
        psi<-unlist(tapply(psi, opz$id.psi.group, sort), use.names =FALSE)
        #browser()
        
        a<-optimize(search.min, c(0,1), psi=psi, psi.old=psi.old, X=XREG, y=y, w=w, offs=offs)
        k.values[length(k.values) + 1] <- use.k <- a$minimum
        L1<- a$objective
        #L1.k[length(L1.k) + 1] <- L1<- a$objective
        psi <- psi*use.k + psi.old* (1-use.k)
        psi<- adj.psi(psi, limZ)
        if (!is.null(digits)) psi <- round(psi, digits)
        PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
        U1 <- (Z - PSI) * (Z > PSI)
        
        
        
        #if (pow[1] != 1) U1 <- U1^pow[1]
        #obj1 <- try(mylm(cbind(XREG, U1), y, w, offs), silent = TRUE)
        #if (class(obj1)[1] == "try-error") obj1 <- try(lm.wfit(cbind(XREG, U1), y, w, offs), silent = TRUE)

        if (visual) {
          flush.console()
          cat(paste("iter = ", sprintf("%2.0f",it),
                    "  dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), L1), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
                    "  k = ", sprintf("%2.3f", use.k),
                    "  n.psi = ",formatC(length(unlist(psi)),digits=0,format="f"), 
                    "  est.psi = ",paste(formatC(unlist(psi),digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                    sep=""), "\n")
        }
        epsilon <- if (conv.psi) 
            max(abs((psi - psi.old)/psi.old))
        else (L0 - L1)/(abs(L0) + 0.1)
        L0 <- L1
        U <- U1
        k.values[length(k.values) + 1] <- use.k
        psi.values[[length(psi.values) + 1]] <- psi
        dev.values[length(dev.values) + 1] <- L0
        id.psi.far <- far.psi(Z, PSI, id.psi.group, TRUE, fc = opz$fc)
        id.psi.in <- in.psi(limZ, PSI, TRUE)
        id.psi.ok <- id.psi.in & id.psi.far
        if (!all(id.psi.ok)) {
            if (fix.npsi) {
                psi <- psi * ifelse(id.psi.far, 1, attr(id.psi.far, 
                  "factor"))
                PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), 
                  ncol = length(psi))
                id.psi.changed[it] <- TRUE
            }
            else {
                Z <- Z[, id.psi.ok, drop = FALSE]
                PSI <- PSI[, id.psi.ok, drop = FALSE]
                rangeZ <- rangeZ[, id.psi.ok, drop = FALSE]
                limZ <- limZ[, id.psi.ok, drop = FALSE]
                nomiOK <- nomiOK[id.psi.ok]
                id.psi.group <- id.psi.group[id.psi.ok]
                psi.old <- psi.old[id.psi.ok]
                psi <- psi[id.psi.ok]
                names(psi) <- id.psi.group
                if (ncol(PSI) <= 0) {
                  warning(paste("All breakpoints have been removed after", 
                    it, "iterations.. returning 0"), call. = FALSE)
                  return(0)
                }
            }
        }
        if (it >= it.max) {
            id.warn <- TRUE
            break
        }
    } #end while..
    ##############################################################################
    if (id.psi.changed[length(id.psi.changed)]) 
        warning(paste("Some psi (", (1:length(psi))[!id.psi.far], 
            ") changed after the last iter.", sep = ""), call. = FALSE)
    if (id.warn) 
        warning(paste("max number of iterations (", it, ") attained", 
            sep = ""), call. = FALSE)
    attr(psi.values, "dev") <- dev.values
    attr(psi.values, "k") <- k.values
    psi <- unlist(tapply(psi, id.psi.group, sort))
    names(psi) <- id.psi.group
    names.coef <- names(obj$coefficients)
    #PSI.old <- PSI
    PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
    #if (sd(PSI - PSI.old) > 0 || id.psi.changed[length(id.psi.changed)]) {
    #browser()
    V <- -(Z > PSI)
    colnames(V) <- paste("V", 1:ncol(V), sep = "")
    
    U <- (Z - PSI) * (Z > PSI)
    for(i in 1:length(RList)){#trasforma le U
      UList[[i]]<- cbind(Zseg[,i], U[, id.psi.group==i])%*%invA.RList[[i]] 
      nomiUList[[i]]<- rep(i, ncol(UList[[i]]) )
    }
    U<-do.call(cbind, UList) #X <- cbind(XREG, U, V)

    colnames(U) <- paste("U", 1:ncol(U), sep = "")
    obj <- try(suppressWarnings(glm.fit(cbind(XREG, U), y = y, offset = offs,
                                        weights = w, family = fam, control = glm.control(maxit = maxit.glm), etastart = eta0)), 
               silent = TRUE)
    L1<- obj$dev    
    #browser()
    obj$coefficients <- c(obj$coefficients, rep(0, ncol(V)))
    #names(obj$coefficients) <- names.coef
    obj$epsilon <- epsilon
    obj$it <- it
    obj <- list(obj = obj, it = it, psi = psi, psi.values = psi.values, X=XREG,
        U = U, V = V, rangeZ = rangeZ, epsilon = epsilon, nomiOK = nomiOK, 
        #SumSquares.no.gap = L1, 
        dev.no.gap=L1,
        id.psi.group = id.psi.group, id.warn = id.warn, 
        constr=list(RList=RList, invAList=invAList, invA.RList=invA.RList, nomiUList =nomiUList))
    #SlopeList <- lapply(1:length(RList), function(i) RList[[i]]%*%beta.c[unlist(nomiUList)==i])
    return(obj)
}
