seg.def.fit<-function (obj, Z, PSI, mfExt, opz, return.all.sol = FALSE) {
    useExp.k = TRUE
    search.min<-function(h, psi, psi.old) { #, mfExt , X, y, w, offs
      psi.ok<- psi*h + psi.old*(1-h)
      PSI <- matrix(rep(psi.ok, rep(n, length(psi.ok))), ncol = length(psi.ok))
      U1 <- (Z - PSI) * (Z > PSI)
      #if (pow[1] != 1) U1 <- U1^pow[1]
      if(is.null(opz$mydesign)){
        for (i in 1:ncol(U1)) {
          mfExt[nomiU[i]] <- U1[, i]
        }
      } else {
        for (i in 1:ncol(U1)) {
          mydesign$variables[nomiU[i]] <- U1[, i]
        }
      }
      obj1 <- suppressWarnings(try(eval(call.noV, envir = mfExt), silent = TRUE))
      L1 <- if (class(obj1)[1] == "try-error") L0 + 10
          else eval(parse(text = fn.obj), list(x = obj1))
      L1
    }
    
    est.k <- function(x1, y1, L0) {
        ax <- log(x1)
        .x <- cbind(1, ax, ax^2)
        b <- drop(solve(crossprod(.x), crossprod(.x, y1)))
        const <- b[1] - L0
        DD <- sqrt(b[2]^2 - 4 * const * b[3])
        kk <- exp((-b[2] + DD)/(2 * b[3]))
        return(round(kk))
    }
    dpmax <- function(x, y, pow = 1) {
        if (pow == 1) 
            -(x > y)
        else -pow * ((x - y) * (x > y))^(pow - 1)
    }
    in.psi <- function(LIM, PSI, ret.id = TRUE) {
        a <- PSI[1, ] < LIM[1, ]
        b <- PSI[1, ] > LIM[2, ]
        is.ok <- !a & !b
        if (ret.id) 
            return(is.ok)
        isOK <- all(is.ok) && all(!is.na(is.ok))
        isOK
    }
    far.psi<-function(Z, PSI, id.psi.group, ret.id=TRUE, fc=.93) {
        #check if psi are far from the boundaries ..s
        #   returns TRUE, if fine.
        #id.far.ok<-sapply(unique(id.psi.group), function(.x) (table(rowSums(((Z>PSI)[,id.psi.group==.x,drop=FALSE])))>=2)[-1]) #[-1] esclude lo zero, x<psi[1] 
        #id.far.ok<-sapply(unique(id.psi.group), function(.x) (tabulate(rowSums(((Z>PSI)[,id.psi.group==.x,drop=FALSE]))+1)>=2)[-1]) #[-1] esclude lo zero, x<psi[1]
        #16/01/20:
        # se un psi assume l'estremo superiore "Z>PSI" non se ne accorge, mentre Z>=PSI, si.. Il contrario e' vero con estremo inf e Z>PSI
        nSeg<-length(unique(id.psi.group))
        npsij<-tapply(id.psi.group,id.psi.group,length)
        nj<-sapply(unique(id.psi.group), function(.x) { tabulate(rowSums((Z>PSI)[,id.psi.group==.x,drop=FALSE])+1) }, simplify = FALSE)    
        ff<-id.far.ok<-vector("list",length=nSeg) 
        for(i in 1:nSeg){
            if(length(nj[[i]])!=npsij[i]+1) nj[[i]]<- tabulate(rowSums((Z>=PSI)[,id.psi.group==i,drop=FALSE])+1)
            id.ok<-(nj[[i]] >= 2)
            id.far.ok[[i]] <- id.ok[-length(id.ok)] & id.ok[-1] #esattamente uguale al numero di psi del gruppo i
            ff[[i]]<-ifelse(diff(nj[[i]])>0, 1/fc, fc)
        }
        id.far.ok<-unlist(id.far.ok)
        ff<-unlist(ff)
        if(!ret.id) {return(all(id.far.ok))
        } else {
            attr(id.far.ok,"factor") <- ff
            return(id.far.ok) 
        }
        #if(ret.id) return(id.far.ok) else return(all(id.far.ok))
    }
    adj.psi <- function(psii, LIM) {
        pmin(pmax(LIM[1, ], psii), LIM[2, ])
    }
    fn.costr <- function(n.psi, isLeft = 1, isInterc = 1) {
        IU <- -diag(n.psi)
        sumU <- diag(n.psi)
        sumU[row(sumU) > col(sumU)] <- 1
        if (isLeft) {
            sumU <- cbind(1, sumU)
            IU <- diag(c(1, -rep(1, n.psi)))
        }
        A <- rbind(IU, sumU)
        if (isInterc) {
            A <- rbind(0, A)
            A <- cbind(c(1, rep(0, nrow(A) - 1)), A)
        }
        A <- cbind(A, matrix(0, nrow(A), n.psi))
        A
    }
    vincoli <- FALSE
    c1 <- apply((Z <= PSI), 2, all)
    c2 <- apply((Z >= PSI), 2, all)
    if (sum(c1 + c2) != 0 || is.na(sum(c1 + c2))) 
        stop("psi out of the range")
    n <- nrow(Z)
    min.step <- opz$min.step
    rangeZ <- apply(Z, 2, range)
    alpha <- opz$alpha
    limZ <- apply(Z, 2, quantile, names = FALSE, probs = c(alpha[1], alpha[2]))
    digits <- opz$digits
    pow <- opz$pow
    nomiOK <- opz$nomiOK
    toll <- opz$toll
    hh <- opz$h
    conv.psi <- opz$conv.psi
    gap <- opz$gap
    stop.if.error <- opz$stop.if.error
    fix.npsi <- opz$fix.npsi
    dev.new <- opz$dev0
    visual <- opz$visual
    id.psi.group <- opz$id.psi.group
    it.max <- old.it.max <- opz$it.max
    fc<-opz$fc
    psi <- PSI[1, ]
    psi<-adj.psi(psi, limZ)
    PSI<-matrix(psi,nrow=n, ncol=ncol(PSI), byrow=TRUE)
    names(psi) <- id.psi.group
    epsilon <- 10
    dev.values <- psi.values <- NULL
    it <- 0
    epsilon <- 10
    k.values <- dev.values <- NULL
    psi.values <- list()
    psi.values[[length(psi.values) + 1]] <- NA
    nomiU <- opz$nomiU
    nomiV <- opz$nomiV
    call.ok <- opz$call.ok
    call.noV <- opz$call.noV
    toll <- opz$toll
    mydesign<- opz$mydesign
    
    #browser()
    
    if (!in.psi(limZ, PSI, FALSE)) 
        stop("starting psi out of the range", call. = FALSE)
    if (!far.psi(Z, PSI, id.psi.group, FALSE)) 
      stop("psi starting values too close each other or at the boundaries. Please change them (e.g. set 'quant=TRUE' 
          in seg.control()), or decrease their number.", call. = FALSE)
    n.psi1 <- ncol(Z)
    if (is.null(opz$constr)) 
        opz$constr <- 0
    if ((opz$constr %in% 1:2) && class(obj)[1] == "rq") {
        vincoli <- TRUE
        call.ok$method <- "fnc"
        call.ok$R <- quote(R)
        call.ok$r <- quote(r)
        call.noV$method <- "fnc"
        call.noV$R <- quote(R.noV)
        call.noV$r <- quote(r)
    }
    fn.obj <- opz$fn.obj
    U <- ((Z - PSI) * (Z > PSI))
    colnames(U) <- nomiU
    if (pow[1] != 1) U <- U^pow[1]
    
   
    
    obj0 <- suppressWarnings(try(eval(call.noV, envir = mfExt), silent = TRUE))
    if ("try-error" %in% class(obj0)) 
        stop("The first fit with U variables does not work..", 
            call. = FALSE)
    L0 <- eval(parse(text = fn.obj), list(x = obj0))
    
    
    if(it.max==0){
      colnames(U) <- paste("U", 1:ncol(U), sep = "")
      V <- -(Z > PSI)
      colnames(V) <- paste("V", 1:ncol(V), sep = "")
      obj <- obj0 #lm.wfit(x = cbind(XREG, U), y = y, w = w, offset = offs)
      L1 <- L0 #sum(obj$residuals^2 * w)
      obj$coefficients <- c(obj$coefficients, rep(0, ncol(V)))
      #names(obj$coefficients) <- names.coef
      obj$epsilon <- epsilon
      obj$it <- it
      obj <- list(obj = obj, it = it, psi = psi, psi.values = psi.values, 
                  U = U, V = V, rangeZ = rangeZ, epsilon = epsilon, nomiOK = nomiOK, 
                  SumSquares.no.gap = L1, id.psi.group = id.psi.group, 
                  id.warn = TRUE, nomiV = nomiV, nomiU = nomiU, mfExt = mfExt, 
                  mydesign=mydesign)
      return(obj)
    }
    

    n.intDev0 <- nchar(strsplit(as.character(L0), "\\.")[[1]][1])
    dev.values[length(dev.values) + 1] <- opz$dev0
    dev.values[length(dev.values) + 1] <- L0
    psi.values[[length(psi.values) + 1]] <- psi
    if (visual) {
        cat(paste("iter = ", sprintf("%2.0f", 0), "  min.f = ", 
            sprintf(paste("%", n.intDev0 + 6, ".5f", sep = ""), 
                L0), "  k = ", sprintf("%2.0f", NA), "  n.psi = ", 
            formatC(length(unlist(psi)), digits = 0, format = "f"), 
            "  ini.psi = ", paste(formatC(unlist(psi), digits = 3, 
                format = "f"), collapse = "  "), sep = ""), "\n")
    }
    id.warn <- FALSE
    id.psi.changed<-rep(FALSE, it.max)
    while (abs(epsilon) > toll) {
        it <- it + 1
        n.psi0 <- n.psi1
        n.psi1 <- ncol(Z)
        if (n.psi1 != n.psi0) {
            U <- ((Z - PSI) * (Z > PSI))
            if (pow[1] != 1) U <- U^pow[1]
            for (i in 1:ncol(U)) {
              mfExt[nomiU[i]] <- U[, i]
            }
            obj0 <- suppressWarnings(try(eval(call.noV, envir = mfExt), silent = TRUE))
            L0 <- eval(parse(text = fn.obj), list(x = obj0))
        }
        V <- dpmax(Z, PSI, pow = pow[2])
        if(is.null(opz$mydesign)){
          for (i in 1:n.psi1) {
            mfExt[nomiU[i]] <- U[, i]
            mfExt[nomiV[i]] <- V[, i]
          }
        } else {
          for (i in 1:n.psi1) {
            mydesign$variables[nomiU[i]] <- U[, i]
            mydesign$variables[nomiV[i]] <- V[, i]
          }
        }
        R <- fn.costr(ncol(U), 1, 1)
        R.noV <- R[, -((ncol(R) - 1) + seq_len(ncol(U))), drop = FALSE]
        r <- rep(0, nrow(R))
        obj <- suppressWarnings(eval(call.ok, envir = mfExt))
        beta.c <- unlist( unique(coef(obj)[nomiU]))
        gamma.c <-unlist(  unique(coef(obj)[nomiV]))
#browser()
        if (any(is.na(c(beta.c, gamma.c)))) {
            if (fix.npsi) {
                if (return.all.sol) 
                  return(list(dev.values, psi.values))
                else stop("breakpoint estimate too close or at the boundary causing NA estimates.. too many breakpoints being estimated?", 
                  call. = FALSE)
            } else {
                id.coef.ok <- !is.na(gamma.c)
                psi <- psi[id.coef.ok]
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
        psi<- adj.psi(psi, limZ)
        
        a<-optimize(search.min, c(0,1), psi=psi, psi.old=psi.old)
        k.values[length(k.values) + 1] <- use.k <- a$minimum
        L1<- a$objective
        psi <- psi*use.k + psi.old* (1-use.k)
        psi<- adj.psi(psi, limZ)
        psi<-unlist(tapply(psi, opz$id.psi.group, sort), use.names =FALSE)
        if (!is.null(digits)) psi <- round(psi, digits)
        #PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
        PSI <- matrix(psi, ncol=length(psi), nrow=n, byrow = TRUE)
        U <- (Z - PSI) * (Z > PSI)
        #if (pow[1] != 1) U <- U^pow[1]

        if(is.null(opz$mydesign)){
          for(i in 1:ncol(U)) mfExt[nomiU[i]] <- U[, i]
        } else {
          for(i in 1:ncol(U))  mydesign$variables[nomiU[i]] <- U[, i]
        }
        

        # obj1 <- suppressWarnings(try(eval(call.noV, envir = mfExt), silent = TRUE))
        # L1 <- if (class(obj1)[1] == "try-error") 
        #     L0 + 10
        # else eval(parse(text = fn.obj), list(x = obj1))
        # use.k <- k <- 1
        # L1.k <- NULL
        # L1.k[length(L1.k) + 1] <- L1
        # # while (L1 > L0) {
        # #     k <- k + 1
        # #     use.k <- if (useExp.k) 
        # #         2^(k - 1)
        # #     else k
        # #     psi <- psi.old + (gamma.c/beta.c)/(use.k * h)
        # #     if (!is.null(digits)) 
        # #         psi <- round(psi, digits)
        # #     PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
        # #     U <- (Z - PSI) * (Z > PSI)
        # #     if (pow[1] != 1) 
        # #         U <- U^pow[1]
        # #     if(is.null(opz$mydesign)){
        # #       for (i in 1:ncol(U)) mfExt[nomiU[i]] <- U[, i]
        # #       #obj1 <- suppressWarnings(try(eval(call.noV, envir = mfExt), silent = TRUE))
        # #     } else {
        # #       for (i in 1:ncol(U)) mydesign$variables[nomiU[i]] <- U[, i]
        # #       #obj1 <- suppressWarnings(try(eval(call.noV), silent = TRUE))
        # #     }
        # #     obj1 <- suppressWarnings(try(eval(call.noV, envir = mfExt), silent = TRUE))
        # #     L1 <- if (class(obj1)[1] == "try-error") 
        # #         L0 + 10
        # #     else eval(parse(text = fn.obj), list(x = obj1))
        # #     L1.k[length(L1.k) + 1] <- L1
        # #     if (1/(use.k * h) < min.step) {
        # #         break
        # #     }
        # # }
        if (visual) {
            flush.console()
            cat(paste("iter = ", sprintf("%2.0f", it), "  min.f = ", 
                sprintf(paste("%", n.intDev0 + 6, ".5f", sep = ""), 
                  L1), "  k = ", sprintf("%2.3f", use.k), "  n.psi = ", 
                formatC(length(unlist(psi)), digits = 0, format = "f"), 
                "  est.psi = ", paste(formatC(unlist(psi), digits = 3, 
                  format = "f"), collapse = "  "), sep = ""), 
                "\n")
        }
        epsilon <- if (conv.psi) 
            max(abs((psi - psi.old)/psi.old))
        else (L0 - L1)/(abs(L0) + 0.1)
        L0 <- L1
        #U<-U1
        k.values[length(k.values) + 1] <- use.k
        psi.values[[length(psi.values) + 1]] <- psi
        dev.values[length(dev.values) + 1] <- L0
        id.psi.far <- far.psi(Z, PSI, id.psi.group, TRUE, fc=opz$fc)
        id.psi.in <- in.psi(limZ, PSI, TRUE)
        id.psi.ok <- id.psi.in & id.psi.far
        if (!all(id.psi.ok)) {
            if (fix.npsi) {
              #psi <- psi * ifelse(id.psi.far, 1, 0.9)
              psi <- psi * ifelse(id.psi.far, 1, attr(id.psi.far, "factor"))  
              PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
              id.psi.changed[it]<-TRUE
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
    } #end while
    ###############################################################################################
    if(id.psi.changed[length(id.psi.changed)]) warning(paste("Some psi (", (1:length(psi))[!id.psi.far],
                                                             ") changed after the last iter.",sep=""), call. = FALSE)
    if (id.warn)  warning(paste("max number of iterations (", it, ") attained", sep = ""), call. = FALSE)
    attr(psi.values, "dev") <- dev.values
    attr(psi.values, "k") <- k.values
    psi <- unlist(tapply(psi, id.psi.group, sort))
    names(psi) <- id.psi.group
    names.coef <- names(coef(obj))
    PSI.old <- PSI
    PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
    #if(sd(PSI-PSI.old)>0 || id.psi.changed[length(id.psi.changed)]){
        U <- (Z - PSI) * (Z > PSI)
        colnames(U) <- paste("U", 1:ncol(U), sep = "")
        V <- -(Z > PSI)
        colnames(V) <- paste("V", 1:ncol(V), sep = "")
        if(is.null(opz$mydesign)){
          for (i in 1:n.psi1) {
            mfExt[nomiU[i]] <- U[, i]
            mfExt[nomiV[i]] <- V[, i]
          }
        } else {
          mydesign$variables[nomiU[i]] <- U[, i]
          mydesign$variables[nomiV[i]] <- V[, i]
        }
        obj <- suppressWarnings(try(eval(call.noV, envir = mfExt), silent = TRUE))
        L1 <- eval(parse(text = fn.obj), list(x = obj))
    #} else {
    #    obj <- obj1
    #}
    nomeCoef <- grep("coef", names(obj), value = TRUE)
    if(length(nomeCoef)==0){
        nomeCoef <- grep("estimate", names(obj), value = TRUE) 
    }
    if(length(nomeCoef)==0) stop("I can't extract the estimated coefficients")
    if(is.list(obj[[nomeCoef]])) {
        obj[[nomeCoef]][[1]] <- c(obj[[nomeCoef]][[1]], rep(0, ncol(V)))
        names(obj[[nomeCoef]][[1]]) <- names.coef[1:length(obj[[nomeCoef]][[1]])]
        } else {
        nomiconV<- c( names(obj[[nomeCoef]]), sub("V", "psi", nomiV)) 
        obj[[nomeCoef]] <- c(obj[[nomeCoef]], rep(0, ncol(V)))
        #se i coeff includono un altro parametro (ed es., la varianza come per censReg), l'ordine deve essere 
        #  rispettato.. mentre "names.coef" 
        #names(obj[[nomeCoef]]) <- names.coef
        names(obj[[nomeCoef]]) <- nomiconV
        }
    
    obj$epsilon <- epsilon
    obj$it <- it
    obj <- list(obj = obj, it = it, psi = psi, psi.values = psi.values, 
        U = U, V = V, rangeZ = rangeZ, epsilon = epsilon, nomiOK = nomiOK, 
        SumSquares.no.gap = L1, id.psi.group = id.psi.group, 
        id.warn = id.warn, nomiV = nomiV, nomiU = nomiU, mfExt = mfExt, mydesign=mydesign)
    if (vincoli) {
        obj$R <- R
        obj$R.noV <- R.noV
        obj$r <- r
    }
    return(obj)
}
