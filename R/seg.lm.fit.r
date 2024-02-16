seg.lm.fit <-function (y, XREG, Z, PSI, w, offs, opz, return.all.sol = FALSE) 
{
    useExp.k = TRUE
    search.min<-function(h, psi, psi.old, X, y, w, offs) {
        psi.ok<- psi*h + psi.old*(1-h)
        #PSI <- matrix(rep(psi.ok, rep(n, length(psi.ok))), ncol = length(psi.ok))
        PSI <- matrix(psi.ok, nrow=n, ncol = length(psi.ok), byrow=TRUE)
        U1 <- (Z - PSI) * (Z > PSI)
        #if (pow[1] != 1) U1 <- U1^pow[1]
        obj1 <- try(mylm(cbind(X, U1), y, w, offs), silent = TRUE)
        if (class(obj1)[1] == "try-error") obj1 <- try(lm.wfit(cbind(X, U1), y, w, offs), silent = TRUE)
        L1 <- if (class(obj1)[1] == "try-error") L0 + 10
        else sum(obj1$residuals^2 * w)
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
    mylm <- function(x, y, w, offs = rep(0, length(y))) {
        x1 <- x * sqrt(w)
        y <- y - offs
        y1 <- y * sqrt(w)
        b <- drop(solve(crossprod(x1), crossprod(x1, y1)))
        fit <- drop(tcrossprod(x, t(b)))
        r <- y - fit
        o <- list(coefficients = b, fitted.values = fit, residuals = r, 
            df.residual = length(y) - length(b))
        o
    }
    mylmADD <- function(invXtX, X, v, Xty, y) {
        vtv <- sum(v^2)
        Xtv <- crossprod(X, v)
        m <- invXtX %*% Xtv
        d <- drop(1/(vtv - t(Xtv) %*% m))
        r <- -d * m
        invF <- invXtX + d * tcrossprod(m)
        newINV <- rbind(cbind(invF, r), c(t(r), d))
        b <- crossprod(newINV, c(Xty, sum(v * y)))
        fit <- tcrossprod(cbind(X, v), t(b))
        r <- y - fit
        o <- list(coefficients = b, fitted.values = fit, residuals = r)
        o
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
    n <- length(y)
    min.step <- opz$min.step
    rangeZ <- apply(Z, 2, range)
    alpha <- opz$alpha
    limZ <- apply(Z, 2, quantile, names = FALSE, probs = c(alpha[1], alpha[2]))

    #browser()
    
    psi <- PSI[1, ]
    psi<-adj.psi(psi, limZ)
    PSI<-matrix(psi,nrow=n, ncol=ncol(PSI), byrow=TRUE)
    id.psi.group <- opz$id.psi.group
    conv.psi <- opz$conv.psi
    hh <- opz$h
    digits <- opz$digits
    pow <- opz$pow
    nomiOK <- opz$nomiOK
    toll <- opz$toll
    gap <- opz$gap
    fix.npsi <- opz$stop.if.error
    dev.new <- opz$dev0
    visual <- opz$visual
    it.max <- old.it.max <- opz$it.max
    fc <- opz$fc
    names(psi) <- id.psi.group
    it <- 0
    epsilon <- 10
    k.values <- dev.values <- NULL
    psi.values <- list()
    psi.values[[length(psi.values) + 1]] <- NA
    sel.col.XREG <- unique(sapply(colnames(XREG), function(x) match(x, 
        colnames(XREG))))
    if (is.numeric(sel.col.XREG)) 
        XREG <- XREG[, sel.col.XREG, drop = FALSE]
    invXtX <- opz$invXtX
    Xty <- opz$Xty
    
    #browser()
    if (!in.psi(limZ, PSI, FALSE)) 
        stop("starting psi out of the range.. see 'alpha' in seg.control.", 
            call. = FALSE)
    if (!far.psi(Z, PSI, id.psi.group, FALSE)) 
        stop("psi starting values too close each other or at the boundaries. Please change them (e.g. set 'quant=TRUE' 
          in seg.control()), or decrease their number.", call. = FALSE)
    n.psi1 <- ncol(Z)
    U <- ((Z - PSI) * (Z > PSI))
    if (pow[1] != 1) U <- U^pow[1]
    
    if(it.max==0){
      colnames(U) <- paste("U", 1:ncol(U), sep = "")
      V <- -(Z > PSI)
      colnames(V) <- paste("V", 1:ncol(V), sep = "")
      obj <- lm.wfit(x = cbind(XREG, U), y = y, w = w, offset = offs)
      L1 <- sum(obj$residuals^2 * w)
      obj$coefficients <- c(obj$coefficients, rep(0, ncol(V)))
      #names(obj$coefficients) <- names.coef
      obj$epsilon <- epsilon
      obj$it <- it
      obj <- list(obj = obj, it = it, psi = psi, psi.values = psi.values, 
                U = U, V = V, rangeZ = rangeZ, epsilon = epsilon, nomiOK = nomiOK, 
                SumSquares.no.gap = L1, id.psi.group = id.psi.group, 
                id.warn = TRUE)
      return(obj)
    }
    
    

    obj0 <- try(mylm(cbind(XREG, U), y, w, offs), silent = TRUE)
    if (class(obj0)[1] == "try-error") obj0 <- lm.wfit(cbind(XREG, U), y, w, offs)
    L0 <- sum(obj0$residuals^2 * w)
    n.intDev0 <- nchar(strsplit(as.character(L0), "\\.")[[1]][1])
    dev.values[length(dev.values) + 1] <- opz$dev0
    dev.values[length(dev.values) + 1] <- L0
    psi.values[[length(psi.values) + 1]] <- psi
    if (visual) {
        cat(paste("iter = ", sprintf("%2.0f", 0), 
                  "  dev = ",  sprintf(paste("%", n.intDev0 + 6, ".5f", sep = ""), L0), 
                  "  k = ", sprintf("%5.0f", NA), 
                  "  n.psi = ", formatC(length(unlist(psi)), digits = 0, format = "f"), 
                  "  ini.psi = ", paste(formatC(unlist(psi), digits = 3, format = "f"), collapse = "  "), 
                  sep = ""), "\n")
    }
    id.warn <- FALSE
    id.psi.changed <- rep(FALSE, it.max)
    #============================================== inizio ciclo
    while (abs(epsilon) > toll) {
        it <- it + 1
        n.psi0 <- n.psi1
        n.psi1 <- ncol(Z)
        if (n.psi1 != n.psi0) {
            U <- ((Z - PSI) * (Z > PSI))
            if (pow[1] != 1) 
                U <- U^pow[1]
            obj0 <- try(mylm(cbind(XREG, U), y, w, offs), silent = TRUE)
            if (class(obj0)[1] == "try-error") 
                obj0 <- lm.wfit(cbind(XREG, U), y, w, offs)
            L0 <- sum(obj0$residuals^2 * w)
        }
        V <- dpmax(Z, PSI, pow = pow[2])
        X <- cbind(XREG, U, V)
        rownames(X) <- NULL
        colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U", 
            1:ncol(U), sep = ""), paste("V", 1:ncol(V), sep = ""))
        obj <- lm.wfit(x = X, y = y, w = w, offset = offs)
        
        #browser()
        
        beta.c <- coef(obj)[paste("U", 1:ncol(U), sep = "")]
        gamma.c <- coef(obj)[paste("V", 1:ncol(V), sep = "")]
        if (any(is.na(c(beta.c, gamma.c)))) {
            if (fix.npsi) {
                if (return.all.sol) 
                  return(list(dev.values, psi.values))
                else stop("breakpoint estimate too close or at the boundary causing NA estimates.. too many breakpoints being estimated?", 
                  call. = FALSE)
            }
            else {
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
        #psi<- adj.psi(psi, rangeZ) 
        psi<- adj.psi(psi, limZ)
        #browser()
        
        psi<-unlist(tapply(psi, opz$id.psi.group, sort), use.names =FALSE)
        
        a<-optimize(search.min, c(0,1), psi=psi, psi.old=psi.old, X=XREG, y=y, w=w, offs=offs)
        k.values[length(k.values) + 1] <- use.k <- a$minimum
        L1<- a$objective
        #L1.k[length(L1.k) + 1] <- L1<- a$objective
        psi <- psi*use.k + psi.old* (1-use.k)
        psi<- adj.psi(psi, limZ)
        if (!is.null(digits)) psi <- round(psi, digits)
        #PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
        PSI <- matrix(psi, nrow=n, ncol = length(psi), byrow=TRUE)
        U1 <- (Z - PSI) * (Z > PSI)
        
        #if(it==3) browser()
        
        #if (pow[1] != 1) U1 <- U1^pow[1]
        #obj1 <- try(mylm(cbind(XREG, U1), y, w, offs), silent = TRUE)
        #if (class(obj1)[1] == "try-error") obj1 <- try(lm.wfit(cbind(XREG, U1), y, w, offs), silent = TRUE)

        if (visual) {
            flush.console()
            cat(paste("iter = ", sprintf("%2.0f", it), "  dev = ", 
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
                #PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
                PSI <- matrix(psi, nrow=n, ncol = length(psi), byrow=TRUE)
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
    #PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
    PSI <- matrix(psi, nrow=n, ncol = length(psi), byrow=TRUE)
    #if (sd(PSI - PSI.old) > 0 || id.psi.changed[length(id.psi.changed)]) {
    U <- (Z - PSI) * (Z > PSI)
    colnames(U) <- paste("U", 1:ncol(U), sep = "")
    V <- -(Z > PSI)
    colnames(V) <- paste("V", 1:ncol(V), sep = "")
    obj <- lm.wfit(x = cbind(XREG, U), y = y, w = w, offset = offs)
    L1 <- sum(obj$residuals^2 * w)
    #}
    #else {
    #    obj <- obj1
    #}
    obj$coefficients <- c(obj$coefficients, rep(0, ncol(V)))
    names(obj$coefficients) <- names.coef
    obj$epsilon <- epsilon
    obj$it <- it
    obj <- list(obj = obj, it = it, psi = psi, psi.values = psi.values, 
        U = U, V = V, rangeZ = rangeZ, epsilon = epsilon, nomiOK = nomiOK, 
        SumSquares.no.gap = L1, id.psi.group = id.psi.group, 
        id.warn = id.warn)
    return(obj)
}
