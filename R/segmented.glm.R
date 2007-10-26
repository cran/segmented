`segmented.glm` <-
function(obj, seg.Z, psi, control=seg.control() , model.frame=TRUE, ...){
    it.max <- control$it.max
    toll <- control$toll
    visual <- control$visual
    last <- control$last
    objframe <- update(obj, model = TRUE)
    .y <- model.response(objframe$model)
    a <- model.matrix(seg.Z, data = obj$data)
    a <- subset(a, select = colnames(a)[-1])
    .n <- nrow(a)
    Z <- list()
    for (i in colnames(a)) Z[[length(Z) + 1]] <- a[, i]
    name.Z <- names(Z) <- colnames(a)
    if (!is.list(Z) || !is.list(psi) || is.null(names(Z)) ||
        is.null(names(psi)))
        stop("Z and psi have to be *named* list")
    nomiZpsi <- match(names(Z), names(psi))
    if (!identical(length(Z), length(psi)) || any(is.na(nomiZpsi)))
        stop("Length or names of Z and psi do not match")
    dd <- match(names(Z), names(psi))
    nome <- names(psi)[dd]
    psi <- psi[nome]
    a <- sapply(psi, length)
    b <- rep(1:length(a), times = a)
    Znew <- list()
    for (i in 1:length(psi)) Znew[[length(Znew) + 1]] <- rep(Z[i],
        a[i])
    Z <- matrix(unlist(Znew), nrow = .n)
    colnames(Z) <- rep(nome, a)
    psi <- unlist(psi)
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(.n, k)), ncol = k)
    nomiZ <- rep(nome, times = a)
    ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))],
        function(xxx) {
            1:xxx
        })))
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ, sep = ".")
    KK <- new.env()
    if(is.null(dim(obj$data))){
    for (i in 1:ncol(objframe$model)) assign(all.vars(formula(obj))[i],
                eval(parse(text=all.vars(formula(obj))[i])), envir = KK)
            } else {
    for (i in 1:length(all.vars(formula(obj)))) assign(all.vars(formula(obj))[i],
                obj$data[,all.vars(formula(obj))[i]], envir = KK)
        }
    if (it.max == 0) {
        U <- pmax((Z - PSI), 0)
        colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
        nomiU <- paste("U", colnames(U), sep = "")
        for (i in 1:ncol(U)) assign(nomiU[i], U[, i], envir = KK)
        Fo <- update.formula(formula(obj), as.formula(paste(".~.+",
            paste(nomiU, collapse = "+"))))
        obj <- update(obj, formula = Fo, data = KK)
        if (model.frame)
            obj$model <- data.frame(as.list(KK))
        obj$psi <- psi
        return(obj)
    }
    XREG <- model.matrix(objframe)
    fam <- family(obj)
    o <- model.offset(objframe$model)
    w <- model.weights(objframe$model)
    if (is.null(w))
        w <- rep(1, .n)
    if (is.null(o))
        o <- rep(0, .n)
    initial <- psi
    it <- 1
    epsilon <- 10
    obj0 <- obj
    list.obj <- list(obj)
    maxit.glm <- control$maxit.glm
    psi.values<-NULL
    rangeZ<-apply(Z,2,range)
    while (abs(epsilon) > toll) {
        eta0 <- obj$linear.predictors
        U <- pmax((Z - PSI), 0)
        V <- ifelse((Z > PSI), -1, 0)
        dev.old <- obj$dev
        X <- cbind(XREG, U, V)
        rownames(X) <- NULL
        if (ncol(V) == 1)
            colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c("U", "V")
        else colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U",
            1:k, sep = ""), paste("V", 1:k, sep = ""))
        obj <- suppressWarnings(glm.fit(x = X, y = .y, offset = o,
            weights = w, family = fam, control = glm.control(maxit = maxit.glm),
            etastart = eta0))
        dev.new <- obj$dev
        if (visual) {
            if (it == 1)
                cat(0, " ", formatC(dev.old, 3, format = "f"),
                  "", "(No breakpoint(s))", "\n")
            spp <- if (it < 10)
                ""
            else NULL
            cat(it, spp, "", formatC(dev.new, 3, format = "f"),
                "\n")
        }
        epsilon <- (dev.new - dev.old)/dev.old
        obj$epsilon <- epsilon
        it <- it + 1
        obj$it <- it
        class(obj) <- c("segmented", class(obj))
        list.obj[[length(list.obj) + ifelse(last == TRUE, 0,
            1)]] <- obj
        if (k == 1) {
            beta.c <- coef(obj)["U"]
            gamma.c <- coef(obj)["V"]
        }
        else {
            beta.c <- coef(obj)[paste("U", 1:k, sep = "")]
            gamma.c <- coef(obj)[paste("V", 1:k, sep = "")]
        }
        if (it > it.max) break
        psi.values[[length(psi.values)+1]]<-psi.old <- psi
        psi <- psi.old + gamma.c/beta.c
        PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
        a <- apply((Z < PSI), 2, all)
        b <- apply((Z > PSI), 2, all)
        if (sum(a + b) != 0 || is.na(sum(a + b)))
            stop("(Some) estimated psi out of its range")
        obj$psi <- psi
    }
    psi.values[[length(psi.values)+1]]<- psi
    id.warn <- FALSE
    if (it > it.max) {
        warning("max number of iterations attained", call. = FALSE)
        id.warn <- TRUE
    }
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
    colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
    colnames(Vxb) <- paste(ripetizioni, nomiZ, sep = ".")
    nomiU <- paste("U", colnames(U), sep = "")
    nomiVxb <- paste("psi", colnames(Vxb), sep = "")
    for (i in 1:ncol(U)) {
        assign(nomiU[i], U[, i], envir = KK)
        assign(nomiVxb[i], Vxb[, i], envir = KK)
    }
    nnomi <- c(nomiU, nomiVxb)
    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+",
        paste(nnomi, collapse = "+"))))
    objF <- update(obj0, formula = Fo, data = KK)
    Cov <- vcov(objF)
    id <- match(paste("psi", colnames(Vxb), sep = ""), names(coef(objF)))
    vv <- if (length(id) == 1)
        Cov[id, id]
    else diag(Cov[id, id])
    psi <- cbind(initial, psi, sqrt(vv))
    rownames(psi) <- colnames(Cov)[id]
    colnames(psi) <- c("Initial", "Est.", "St.Err")
    objF$rangeZ<-rangeZ
    objF$psi.history<-psi.values
    objF$psi <- psi
    objF$it <- (it - 1)
    objF$epsilon <- epsilon
    objF$call <- match.call()
    objF$nameUV <- list(U = paste("U", colnames(Vxb), sep = ""),
        V = rownames(psi), Z = name.Z)
    objF$id.warn <- id.warn
    if (model.frame)
        objF$mframe <- data.frame(as.list(KK))
    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last)
        list.obj <- list.obj[[length(list.obj)]]
    return(list.obj)
}
