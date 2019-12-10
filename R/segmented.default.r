segmented.default<-function (obj, seg.Z, psi, npsi, control = seg.control(), model = TRUE, 
    keep.class = FALSE, ...) {

#if("|" %in% all.names(formula(obj))) {
#    nomeY<-all.vars(formula(obj))[1] #nome Y
#    nomiX<-strsplit(as.character(formula(obj))[3],"\\|")[[1]][1]
#    nomiX.disp <-strsplit(as.character(formula(obj))[3],"\\|")[[1]][2]
#    Fo.charac <- paste(nomeY,nomiX,sep="~")
#    Fo <- as.formula(Fo.charac)
#    Fo.conDisp <- as.formula(paste(Fo.charac,nomiX.disp,sep="|"))
    #}

    update.formula1<-function(old,new,...,opt=1){
    #se old e' una formula che contiene "|", questa funzione aggiorna old con new, 
    #       se opt=1 "new" viene inclusa solo nella prima parte e la formula restuita contiene "|"
    #       se opt=2, la parte dopo |, viene aggiunta insieme a "new" e quindi la formula restituita NON contiene | 
        if("|" %in% all.names(old)) {
            nomeY<-all.vars(old)[1] #nome Y
            nomiX<-strsplit(as.character(old)[3],"\\|")[[1]][1]
            nomiX.disp <-strsplit(as.character(old)[3],"\\|")[[1]][2]
            if(opt==2){
                nomiX.all<-paste(nomiX, nomiX.disp,sep="+")
                Fo.charac <- paste(nomeY,nomiX.all,sep="~")
                Fo <- as.formula(Fo.charac)
                    } else {
                Fo.charac <- paste(nomeY,nomiX,sep="~")
                Fo <- as.formula(paste(Fo.charac,nomiX.disp,sep="|"))
                }
            return(Fo)
        } else {
            update.formula(old,new,...)
        }
    }

    dpmax <- function(x, y, pow = 1) {
        if (pow == 1) 
            -(x > y)
        else -pow * ((x - y) * (x > y))^(pow - 1)
    }
    if (is.null(control$fn.obj)) 
        fn.obj <- "-as.numeric(logLik(x))"
    else fn.obj <- control$fn.obj
    if (missing(seg.Z)) {
        if (length(all.vars(formula(obj))) == 2) 
            seg.Z <- as.formula(paste("~", all.vars(formula(obj))[2]))
        else stop("please specify 'seg.Z'")
    }
    n.Seg <- length(all.vars(seg.Z))
    id.npsi <- FALSE
    if (missing(psi)) {
        if (n.Seg == 1) {
            if (missing(npsi)) 
                npsi <- 1
            npsi <- lapply(npsi, function(.x) .x)
            if (length(npsi) != length(all.vars(seg.Z))) 
                stop("seg.Z and npsi do not match")
            names(npsi) <- all.vars(seg.Z)
        }
        else {
            if (missing(npsi)) 
                stop(" with multiple segmented variables in seg.Z, 'psi' or 'npsi' should be supplied", 
                  call. = FALSE)
            if (length(npsi) != n.Seg) 
                stop(" 'npsi' and seg.Z should have the same length")
            if (!all(names(npsi) %in% all.vars(seg.Z))) 
                stop(" names in 'npsi' and 'seg.Z' do not match")
        }
        psi <- lapply(npsi, function(.x) rep(NA, .x))
        id.npsi <- TRUE
    }
    else {
        if (n.Seg == 1) {
            if (!is.list(psi)) {
                psi <- list(psi)
                names(psi) <- all.vars(seg.Z)
            }
        }
        else {
            if (!is.list(psi)) 
                stop("with multiple terms in `seg.Z', `psi' should be a named list")
            if (n.Seg != length(psi)) 
                stop("A wrong number of terms in `seg.Z' or `psi'")
            if (!all(names(psi) %in% all.vars(seg.Z))) 
                stop("Names in `seg.Z' and `psi' do not match")
        }
    }
    min.step <- control$min.step
    alpha <- control$alpha
    it.max <- old.it.max <- control$it.max
    digits <- control$digits
    toll <- control$toll
    if (toll < 0) 
        stop("Negative tolerance ('tol' in seg.control()) is meaningless", 
            call. = FALSE)
    visual <- control$visual
    stop.if.error <- fix.npsi <- control$fix.npsi
    n.boot <- control$n.boot
    size.boot <- control$size.boot
    gap <- control$gap
    random <- control$random
    pow <- control$pow
    conv.psi <- control$conv.psi
    visualBoot <- FALSE
    if (n.boot > 0) {
        if (!is.null(control$seed)) {
            set.seed(control$seed)
            employed.Random.seed <- control$seed
        }
        else {
            employed.Random.seed <- eval(parse(text = paste(sample(0:9, 
                size = 6), collapse = "")))
            set.seed(employed.Random.seed)
        }
        if (visual) {
            visual <- FALSE
            visualBoot <- TRUE
        }
        if (!stop.if.error) 
            stop("Bootstrap restart only with a fixed number of breakpoints")
    }
    last <- control$last
    K <- control$K
    h <- control$h
    orig.call <- Call <- mf <- obj$call
    orig.call$formula <- mf$formula <- formula(obj)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if (class(mf$formula) == "name" && !"~" %in% paste(mf$formula)) 
        mf$formula <- eval(mf$formula)
    #mf$formula <- update.formula(mf$formula, paste(seg.Z, collapse = ".+"))
    mf$formula <- update.formula1(mf$formula, paste(seg.Z, collapse = ".+"), opt=2)
    mfExt <- mf
    
    
    
    if (!is.null(obj$call$offset) || !is.null(obj$call$weights) || 
        !is.null(obj$call$subset) || !is.null(obj$call$id)) {
        mfExt$formula <- update.formula(mf$formula, paste(".~.+", 
            paste(c(all.vars(obj$call$offset), all.vars(obj$call$weights), 
                all.vars(obj$call$subset), all.vars(obj$call$id)), 
                collapse = "+")))
    }
    mf <- eval(mf, parent.frame())
    n <- nrow(mf)
    nomiOff <- setdiff(all.vars(formula(obj)), names(mf))
    if (length(nomiOff) >= 1) 
        mfExt$formula <- update.formula(mfExt$formula, paste(".~.+", 
            paste(nomiOff, collapse = "+"), sep = ""))
    nomiTUTTI <- all.vars(mfExt$formula)
    nomiNO <- NULL
    for (i in nomiTUTTI) {
        r <- try(eval(parse(text = i), parent.frame()), silent = TRUE)
        if (class(r) != "try-error" && length(r) == 1 && !is.function(r)) 
            nomiNO[[length(nomiNO) + 1]] <- i
        }
    if (!is.null(nomiNO)) 
        mfExt$formula <- update.formula(mfExt$formula, paste(".~.-", 
            paste(nomiNO, collapse = "-"), sep = ""))
    mfExt <- eval(mfExt, parent.frame())
    if (inherits(obj, "coxph")) {
        is.Surv <- NA
        rm(is.Surv)
        for (i in 1:ncol(mfExt)) {
            if (is.Surv(mfExt[, i])) 
                aa <- mfExt[, i][, 1:ncol(mfExt[, i])]
        }
        mfExt <- cbind(aa, mfExt)
    }
    id.seg <- match(all.vars(seg.Z), names(mfExt))
    name.Z <- names(mfExt)[id.seg]
    Z <- mfExt[, id.seg, drop = FALSE]
    n.psi <- length(unlist(psi))
    if (ncol(Z) == 1 && is.vector(psi) && (is.numeric(psi) || 
        is.na(psi))) {
        psi <- list(as.numeric(psi))
        names(psi) <- name.Z
    }
    if (!is.list(psi) || is.null(names(psi))) 
        stop("psi should be a *named* list")
    id.nomiZpsi <- match(colnames(Z), names(psi))
    if ((ncol(Z) != length(psi)) || any(is.na(id.nomiZpsi))) 
        stop("Length or names of Z and psi do not match")
    nome <- names(psi)[id.nomiZpsi]
    psi <- psi[nome]
    if (id.npsi) {
        for (i in 1:length(psi)) {
            K <- length(psi[[i]])
            if (any(is.na(psi[[i]]))) 
                psi[[i]] <- if (control$quant) {
                  quantile(Z[, i], prob = seq(0, 1, l = K + 2)[-c(1, 
                    K + 2)], names = FALSE)
                }
                else {
                  (min(Z[, i]) + diff(range(Z[, i])) * (1:K)/(K + 
                    1))
                }
        }
    }
    else {
        for (i in 1:length(psi)) {
            if (any(is.na(psi[[i]]))) 
                psi[[i]] <- if (control$quant) {
                  quantile(Z[, i], prob = seq(0, 1, l = K + 2)[-c(1, 
                    K + 2)], names = FALSE)
                }
                else {
                  (min(Z[, i]) + diff(range(Z[, i])) * (1:K)/(K + 
                    1))
                }
        }
    }
    initial.psi <- psi
    a <- sapply(psi, length)
    id.psi.group <- rep(1:length(a), times = a)
    Z <- matrix(unlist(mapply(function(x, y) rep(x, y), Z, a, 
        SIMPLIFY = TRUE)), nrow = n)
    colnames(Z) <- nomiZ.vett <- rep(nome, times = a)
    psi <- unlist(psi)
    psi <- unlist(tapply(psi, id.psi.group, sort))
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
    c1 <- apply((Z <= PSI), 2, all)
    c2 <- apply((Z >= PSI), 2, all)
    if (sum(c1 + c2) != 0 || is.na(sum(c1 + c2))) 
        stop("starting psi out of the admissible range")
    ripetizioni <- as.vector(unlist(tapply(id.psi.group, id.psi.group, 
        function(x) 1:length(x))))
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ.vett, sep = ".")
    nomiV <- paste("V", ripetizioni, sep = "")
    nomiV <- paste(nomiV, nomiZ.vett, sep = ".")
    nnomi <- c(nomiU, nomiV)
    U <- (Z - PSI) * (Z > PSI)
    if (pow[1] != 1) 
        U <- U^pow[1]
    colnames(U) <- nomiU
    V <- -(Z > PSI)
    for (i in 1:k) {
        mfExt[nomiU[i]] <- U[, i]
        mfExt[nomiV[i]] <- V[, i]
    }
    Fo <- update.formula1(formula(obj), as.formula(paste(".~.+", 
        paste(nnomi, collapse = "+"))), opt=1)
    Fo.noV <- update.formula1(formula(obj), as.formula(paste(".~.+", 
        paste(nomiU, collapse = "+"))), opt=1)
    call.ok <- update(obj, Fo, evaluate = FALSE, data = mfExt)
    call.noV <- update(obj, Fo.noV, evaluate = FALSE, data = mfExt)
    if (it.max == 0) {
        if (!is.null(call.noV[["subset"]])) 
            call.noV[["subset"]] <- NULL
        obj1 <- eval(call.noV, envir = mfExt)
        return(obj1)
    }
    initial <- psi
    obj0 <- obj
    dev0 <- eval(parse(text = fn.obj), list(x = obj))
    if (length(dev0) <= 0) 
        stop("error in the objective to be minimized, see 'fn.obj' in ?seg.control")
    if (length(dev0) > 1) 
        stop("the objective to be minimized is not scalar, see 'fn.obj' in ?seg.control")
    if (is.na(dev0)) 
        dev0 <- 10
    list.obj <- list(obj)
    nomiOK <- nomiU
    opz <- list(toll = toll, h = h, stop.if.error = stop.if.error, 
        dev0 = dev0, visual = visual, it.max = it.max, nomiOK = nomiOK, 
        id.psi.group = id.psi.group, gap = gap, visualBoot = visualBoot, 
        pow = pow, digits = digits, conv.psi = conv.psi, alpha = alpha, 
        fix.npsi = fix.npsi, min.step = min.step)
    opz$call.ok <- call.ok
    opz$call.noV <- call.noV
    opz$formula.orig <- formula(obj)
    opz$nomiU <- nomiU
    opz$nomiV <- nomiV
    opz$fn.obj <- fn.obj
    opz <- c(opz, ...)
    if (n.boot <= 0) {
        obj <- seg.def.fit(obj, Z, PSI, mfExt, opz)
    }
    else {
        obj <- seg.def.fit.boot(obj, Z, PSI, mfExt, opz, n.boot = n.boot, 
            size.boot = size.boot, random = random)
    }
    if (!is.list(obj)) {
        warning("No breakpoint estimated", call. = FALSE)
        return(obj0)
    }
    if (!is.null(obj$obj$df.residual) && !is.na(obj$obj$df.residual)) {
        if (obj$obj$df.residual == 0) 
            warning("no residual degrees of freedom (other warnings expected)", 
                call. = FALSE)
    }
    id.psi.group <- obj$id.psi.group
    nomiU <- nomiOK <- obj$nomiOK
    nomiVxb <- sub("U", "psi", nomiOK)
    nomiFINALI <- unique(sub("U[1-9]*[0-9].", "", nomiOK))
    nomiSenzaPSI <- setdiff(name.Z, nomiFINALI)
    if (length(nomiSenzaPSI) >= 1) 
        warning("no breakpoints found for: ", paste(nomiSenzaPSI, 
            " "), call. = FALSE)
    it <- obj$it
    psi <- obj$psi
    psi.values <- if (n.boot <= 0) 
        obj$psi.values
    else obj$boot.restart
    
    U <- obj$U
    V <- obj$V
    id.warn <- obj$id.warn
    for (jj in colnames(V)) {
        VV <- V[, which(colnames(V) == jj), drop = FALSE]
        sumV <- abs(rowSums(VV))
        if (any(table(sumV) <= 1) && stop.if.error) 
            stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
    }
    rangeZ <- obj$rangeZ
    mfExt <- obj$mfExt
    names(mfExt)[match(obj$nomiV, names(mfExt))] <- nomiVxb
    R <- obj$R
    R.noV <- obj$R.noV
    r <- obj$r
    obj <- obj$obj
    k <- length(psi)
    
    #browser()
    
    #coef(obj) ha gia i nomi corretti... 
    #all.coef <- coef(obj)
    #names(all.coef) <- c(names(coef(obj0)), nomiU, nomiVxb)
    #beta.c <- all.coef[nomiU]
    beta.c<-coef(obj)[nomiU]
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
    nnomi <- c(nomiU, nomiVxb)
    for (i in 1:ncol(U)) {
        mfExt[nomiU[i]] <- mf[nomiU[i]] <- U[, i]
        mfExt[nomiVxb[i]] <- mf[nomiVxb[i]] <- Vxb[, i]
    }
    Fo <- update.formula1(formula(obj0), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))), opt=1)
    objF <- update(obj0, Fo, evaluate = FALSE, data = mfExt)
    if (!is.null(objF[["subset"]])) 
        objF[["subset"]] <- NULL
    if (is.null(opz$constr)) 
        opz$constr <- 0
    if ((opz$constr %in% 1:2) && class(obj0) == "rq") {
        objF$method <- "fnc"
        objF$R <- quote(R)
        objF$r <- quote(r)
    }
    objF <- eval(objF, envir = mfExt)
    objF$offset <- obj0$offset
    isNAcoef <- any(is.na(coef(objF)))
    if (isNAcoef) {
        if (stop.if.error) {
            cat("breakpoint estimate(s):", as.vector(psi), "\n")
            stop("at least one coef is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)", 
                call. = FALSE)
        }
        else {
            warning("some estimate is NA: too many breakpoints? 'var(hat.psi)' cannot be computed \n ..returning a 'lm' model", 
                call. = FALSE)
            Fo <- update.formula1(formula(obj0), as.formula(paste(".~.+", 
                paste(nomiU, collapse = "+"))), opt=1)
            objF <- if ((opz$constr %in% 1:2) && class(obj0) == 
                "rq") {
                update(obj0, formula = Fo, R = R.noV, r = r, 
                  method = "fnc", evaluate = TRUE, data = mfExt)
            }
            else {
                update(obj0, Fo, evaluate = TRUE, data = mfExt)
            }
            names(psi) <- nomiVxb
            objF$psi <- psi
            return(objF)
        }
    }
    
    #4/12/19: modifica fatta per consentire betareg.. Attenzione
    #semplicemente controlla se la componente "coef*" e' una lista o no..
    #COSA succede con geese models?
    if(!is.list(objF[[grep("coef", names(objF), value = TRUE)]])){
            names.coef <- names(coef(objF))
            names(obj[[grep("coef", names(obj), value = TRUE)]]) <- names(objF[[grep("coef", names(objF), value = TRUE)]])
            objF[[grep("coef", names(objF), value = TRUE)]][names.coef] <- coef(obj)[names.coef]
        } else {
            #names.coef <- names(objF[[grep("coef", names(objF), value = TRUE)]][[1]])
            names(obj[[grep("coef", names(obj), value = TRUE)]][[1]]) <- names(objF[[grep("coef", names(objF), value = TRUE)]][[1]])
            objF[[grep("coef", names(objF), value = TRUE)]][[1]] <- 
                obj[[grep("coef", names(obj), value = TRUE)]][[1]]
            objF[[grep("coef", names(objF), value = TRUE)]][[2]] <- 
                obj[[grep("coef", names(obj), value = TRUE)]][[2]]
        }

    if (!is.null(objF$pseudo.r.squared)) 
        objF$pseudo.r.squared <- obj$pseudo.r.squared
    if (!is.null(objF$geese$beta)) 
        objF$geese$beta <- obj$coefficients #oppure objF$coefficients?
    if (!is.null(objF$geese$gamma)) 
        objF$geese$gamma <- obj$geese$gamma
    if (!is.null(objF$geese$alpha)) 
        objF$geese$alpha <- obj$geese$alpha
    if (!is.null(objF$fitted.values)) 
        objF$fitted.values <- obj$fitted.values
    if (!is.null(objF$residuals)) 
        objF$residuals <- obj$residuals
    if (!is.null(objF$linear.predictors)) 
        objF$linear.predictors <- obj$linear.predictors
    if (!is.null(objF$deviance)) 
        objF$deviance <- obj$deviance
    if (!is.null(objF$weights)) 
        objF$weights <- obj$weights
    if (!is.null(objF$aic)) 
        objF$aic <- obj$aic + 2 * k
    if (!is.null(objF$loglik)) 
        objF$loglik <- obj$loglik
    if (!is.null(objF$rho)) 
        objF$rho <- obj$rho
    if (!is.null(objF$dual)) 
        objF$dual <- obj$dual
    if (!is.null(objF$penalized.deviance)) 
        objF$penalized.deviance <- obj$penalized.deviance
    if (!is.null(objF$ModifiedScores)) 
        objF$ModifiedScores <- c(obj$ModifiedScores, rep(0, k))
    Cov <- try(vcov(objF), silent = TRUE)
    if(inherits(Cov, "try-error")){
    #if (class(Cov) == "try-error") { 
        warning("cannot compute the covariance matrix", call. = FALSE)
        vv <- NA
    }
    else {
        vv <- Cov[nomiVxb, nomiVxb, drop=FALSE]
    }
    ris.psi <- matrix(NA, length(psi), 3)
    colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
    rownames(ris.psi) <- nomiVxb
    ris.psi[, 2] <- psi
    ris.psi[, 3] <- sqrt(diag(vv))
    a <- tapply(id.psi.group, id.psi.group, length)
    a.ok <- NULL
    for (j in name.Z) {
        if (j %in% nomiFINALI) {
            a.ok[length(a.ok) + 1] <- a[1]
            a <- a[-1]
        }
        else {
            a.ok[length(a.ok) + 1] <- 0
        }
    }
    initial <- unlist(mapply(function(x, y) {
        if (is.na(x)[1]) 
            rep(x, y)
        else x
    }, initial.psi[nomiFINALI], a.ok[a.ok != 0], SIMPLIFY = TRUE))
    if (opz$stop.if.error) 
        ris.psi[, 1] <- initial
    objF$rangeZ <- rangeZ
    objF$psi.history <- psi.values
    objF$psi <- ris.psi
    objF$it <- it
    objF$epsilon <- obj$epsilon
    objF$call <- match.call()
    objF$nameUV <- list(U = drop(nomiU), V = rownames(ris.psi), 
        Z = nomiFINALI)
    objF$id.group <- if (length(name.Z) <= 1) 
        -rowSums(as.matrix(V))
    objF$id.psi.group <- id.psi.group
    objF$id.warn <- id.warn
    objF$orig.call <- orig.call
    if (model) 
        objF$model <- mf
    if (n.boot > 0) 
        objF$seed <- employed.Random.seed
    if (keep.class) 
        class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last) 
        list.obj <- list.obj[[length(list.obj)]]
    warning("The returned fit is ok, but not of class 'segmented'. If interested, call explicitly the segmented methods (plot.segmented, confint.segmented,..)", 
        call. = FALSE)
    return(list.obj)
}
