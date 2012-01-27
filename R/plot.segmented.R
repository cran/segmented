plot.segmented<-function (x, term = NULL, se = FALSE, const = coef(x)["(Intercept)"],
    add = FALSE, linkinv = FALSE, show.gap = TRUE, rev.sgn=FALSE, n.points=10, ...){
    if(is.na(const)) stop("Undefined `const' argument.. a fitted model without intercept?")
    if (se) {
        se <- FALSE
        warning("se=TRUE not (yet) implemented", call. = FALSE)
    }
    colori<-list(...)$col
    if(length(colori)<=0) colori<-1
    lwds<-list(...)$lwd
    if(length(lwds)<=0) lwds<-1
    ltys<-list(...)$lty
    if(length(ltys)<=0) ltys<-1
    
    mio.plot <- function(..., col, lty, lwd) plot.default(...)
    mio.lines <- function(..., col, ylim, xlim, pch, xlab, ylab) lines(..., col=col)
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
    id.term <- c(match(term,names(coef(x))), 
        grep(paste("\\.",term,"$",sep=""), names(coef(x)), value = FALSE))
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
    idpsi <- grep(paste("\\.",term,"$",sep=""), rownames(x$psi), value = FALSE)
    psi <- x$psi[idpsi, 2]
    xx <- seq(min(x$rangeZ[, term]), max(x$rangeZ[, term]), length = n.points)
    xx<- sort(c(xx,psi,(psi+psi/10000)))
    #xx <- sort(c(xx,psi,psi+diff(range(xx))/300,psi-diff(range(xx))/300))
    #xx <- sort(c(xx, psi, psi + .Machine$double.eps*diff(range(xx)),
    #          psi - .Machine$double.neg.eps*diff(range(xx))))
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
    } else {        
    yhat <- yhat + const
    }
    ylab <- "Fitted Values (on the link scale)"
    if (inherits(x, what = "glm", which = FALSE) && linkinv) {
        yhat <- x$family$linkinv(yhat)
        if (se) {
            inf.yhat <- x$family$linkinv(inf.yhat)
            sup.yhat <- x$family$linkinv(sup.yhat)
        }
        ylab <- "Fitted Values"
    }
    sel.g <- rowSums(V)
    sel.g<-abs(sel.g)+1
    if(length(colori)<length(unique(sel.g))) colori<-rep(colori, length.out=length(unique(sel.g)))
    if(length(lwds)<length(unique(sel.g))) lwds<-rep(lwds, length.out=length(unique(sel.g)))
    if(length(ltys)<length(unique(sel.g))) ltys<-rep(ltys, length.out=length(unique(sel.g)))
    
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
    #for (i in sel.g) mio.lines(xx[sel.g == i], yhat[sel.g ==i], ...)
    for (i in sel.g) {
        list.lines<-c(list(x=xx[sel.g == i], y=yhat[sel.g ==i] ) , list(...))
        list.lines$col<-colori[i]
        list.lines$lwd<-lwds[i]
        list.lines$lty<-ltys[i]
        do.call(lines, list.lines)
        }
    
    if (se) {
        for (i in sel.g) mio.lines(xx[sel.g == i], inf.yhat[sel.g ==
            i], ...)
        for (i in sel.g) mio.lines(xx[sel.g == i], sup.yhat[sel.g ==
            i], ...)
    }
}


