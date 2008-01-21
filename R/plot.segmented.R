plot.segmented<-function (x, term = NULL, se = FALSE, const = coef(x)["(Intercept)"],
    add = FALSE, linkinv = FALSE, show.gap = TRUE, rev.sgn=FALSE, n.points=10, ...){
    if(is.na(const)) stop("Undefined `const' argument")
    if (se) {
        se <- FALSE
        warning("se=TRUE not (yet) implemented", call. = FALSE)
    }
    mio.plot <- function(..., col, lty, lwd) plot.default(...)
    mio.lines <- function(..., ylim, xlim, pch, xlab, ylab) lines(...)
    if (missing(term)) {
        if (length(x$nameUV$Z) > 1) {
            stop("please, specify `term'")
        }
        else {
            term <- x$nameUV$Z
        }
    } else {
      if(!term%in%x$nameUV$Z) stop("invalid `term'")
      }
    id.Z <- match(x$nameUV$Z, names(coef(x)))
    id.U <- match(x$nameUV$U, names(coef(x)))
    id.V <- match(x$nameUV$V, names(coef(x)))
    id.term <- grep(term, names(coef(x)), extended = FALSE)
    id.U <- id.term[id.term %in% id.U]
    id.V <- id.term[id.term %in% id.V]
    id.Z <- id.term[id.term %in% id.Z]
    cof <- c(coef(x)[id.U], coef(x)[id.V] * coef(x)[id.U])
    cof <- if (is.na(id.Z) || length(id.Z) == 0)
        c(0, cof)
    else c(coef(x)[id.Z], cof)
    if (se) {
        A <- diag(c(1, 1, coef(x)[id.U]))
        id.ok <- c(id.Z, id.U, id.V)
        cof <- A %*% coef(x)[id.ok]
        var.cof <- A %*% vcov(x)[id.ok, id.ok] %*% A
    }
    idpsi <- grep(term, rownames(x$psi), extended = FALSE)
    psi <- x$psi[idpsi, 2]
    xx <- seq(min(x$rangeZ[, term]), max(x$rangeZ[, term]), length = n.points)
    #xx <- sort(c(xx,psi,psi+diff(range(xx))/300,psi-diff(range(xx))/300))
    xx <- sort(c(xx, psi, psi + .Machine$double.eps*diff(range(xx)),
              psi - .Machine$double.neg.eps*diff(range(xx))))
    #double.neg.eps
    #if(rev.sgn) xx<- sort(-xx)
    nn <- length(xx)
    k <- length(psi)
    Z <- matrix(rep(xx, times = k), nrow = nn, byrow = FALSE)
    PSI <- matrix(rep(psi, rep(nn, k)), ncol = k)
    U <- pmax((Z - PSI), 0)
    V <- ifelse(Z > PSI, -1, 0)
    X <- cbind(xx, U, V)
    yhat <- drop(X %*% cof)
    if (!show.gap) {
        X <- cbind(xx, U)
        cof <- lm.fit(x = X, y = yhat)$coefficients
        yhat <- drop(X %*% cof)
    }
    if (se) {
        se.yhat <- sqrt(diag(X %*% tcrossprod(var.cof, X)))
        inf.yhat <- yhat - 1.96 * se.yhat
        sup.yhat <- yhat + 1.96 * se.yhat
    }
    else {
        yhat <- yhat + const
    }
    ylab <- "link(Fitted Values)"
    if (inherits(x, what = "glm", which = FALSE) && linkinv) {
        yhat <- x$family$linkinv(yhat)
        if (se) {
            inf.yhat <- x$family$linkinv(inf.yhat)
            sup.yhat <- x$family$linkinv(sup.yhat)
        }
        ylab <- "Fitted Values"
    }
    sel.g <- rowSums(V)
    if(rev.sgn) xx<- -xx
    if (!add) {
        if (!se) {
            mio.plot(xx, yhat, type = "n", ylab = ylab, xlab = term,
                ...)
        }
        else {
            matplot(xx, cbind(yhat, inf.yhat, sup.yhat), type = "n",
                ylab = ylab, xlab = term, ...)
        }
    }
    for (i in sel.g) mio.lines(xx[sel.g == i], yhat[sel.g ==
        i], ...)
    if (se) {
        for (i in sel.g) mio.lines(xx[sel.g == i], inf.yhat[sel.g ==
            i], ...)
        for (i in sel.g) mio.lines(xx[sel.g == i], sup.yhat[sel.g ==
            i], ...)
    }
}


