plot.segmented<-function (x, term, add = FALSE, res = FALSE, conf.level = 0,
    link = TRUE, res.col = 1, rev.sgn = FALSE, const = 0, shade=FALSE, rug=TRUE,
    show.gap=FALSE, ...){
#funzione plot.segmented che consente di disegnare anche i pointwise CI
#Eliminato l'uso di broken.line()
#Basata sulla funzione predict.segmented()
#Ultimo aggiornamento: 5/11/2013
    if(conf.level< 0 || conf.level>.9999) stop("meaningless 'conf.level'")
    show.gap<-FALSE
    if (missing(term)) {
        if (length(x$nameUV$Z) > 1) {
            stop("please, specify `term'")
        }
        else {
            term <- x$nameUV$Z
        }
    }
    else {
        if (!term %in% x$nameUV$Z)
            stop("invalid `term'")
    }
    linkinv <- !link
    opz <- list(...)
    cols <- opz$col
    if (length(cols) <= 0)
        cols <- 1
    lwds <- opz$lwd
    if (length(lwds) <= 0)
        lwds <- 1
    ltys <- opz$lty
    if (length(ltys) <= 0)
        ltys <- 1
    cexs <- opz$cex
    if (length(cexs) <= 0)
        cexs <- 1
    pchs <- opz$pch
    if (length(pchs) <= 0)
        pchs <- 1
    xlabs <- opz$xlab
    if (length(xlabs) <= 0)
        xlabs <- term
    ylabs <- opz$ylab
    if (length(ylabs) <= 0)
        ylabs <- paste("Effect  of ", term, sep = " ")
    a <- intercept(x, term, gap = show.gap)[[1]][, "Est."]
    b <- slope(x, term)[[1]][, "Est."]
    id <- grep(paste("\\.", term, "$", sep = ""), rownames(x$psi),
        value = FALSE)
    est.psi <- x$psi[id, "Est."]
    K <- length(est.psi)
    val <- sort(c(est.psi, x$rangeZ[, term]))
    #---------aggiunta per gli IC
    rangeCI<-NULL
    tipo<-"response"
    if(inherits(x, what = "glm", which = FALSE) && link) tipo<-"link"
    if(conf.level>0) {
        vall<-sort(c(seq(min(val), max(val), l=120), est.psi))
        ciValues<-predict.segmented(x, newdata=vall, se.fit=TRUE, type=tipo, level=conf.level)
        k.alpha<-if(inherits(x, what = "glm", which = FALSE)) abs(qnorm((1-conf.level)/2)) else abs(qt((1-conf.level)/2,ciValues$df))
        ciValues<-cbind(ciValues$fit, ciValues$fit- k.alpha*ciValues$se.fit, ciValues$fit + k.alpha*ciValues$se.fit)
        rangeCI<-range(ciValues)
        #ciValues  è una matrice di length(val)x3. Le 3 colonne: stime, inf, sup
        #polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])), col = "gray", border=NA)
        }

    #---------
    a.ok <- c(a[1], a)
    b.ok <- c(b[1], b)
    y.val <- a.ok + b.ok * val + const
    a.ok1 <- c(a, a[length(a)])
    b.ok1 <- c(b, b[length(b)])
    y.val <- y.val1 <- a.ok1 + b.ok1 * val + const
    s <- 1:(length(val) - 1)
    xvalues <- x$model[, term]
    if (rev.sgn) {
        val <- -val
        xvalues <- -xvalues
    }
    m <- cbind(val[s], y.val1[s], val[s + 1], y.val[s + 1])
    if (inherits(x, what = "glm", which = FALSE) && linkinv) {
        fit <- if (res)
            predict.segmented(x, ifelse(rep(rev.sgn, length(xvalues)),-xvalues,xvalues), type=tipo) + resid(x, "response") + const
            #broken.line(x, term, gap = show.gap, link = link) + resid(x, "response") + const
        else x$family$linkinv(c(y.val, y.val1))
        xout <- sort(c(seq(val[1], val[length(val)], l = 120), val[-c(1, length(val))]))
        l <- approx(as.vector(m[, c(1, 3)]), as.vector(m[, c(2, 4)]), xout = xout)
        id.group <- cut(l$x, val, FALSE, TRUE)
        yhat <- l$y
        xhat <- l$x
        m[, c(2, 4)] <- x$family$linkinv(m[, c(2, 4)])
        if (!add) {
            plot(as.vector(m[, c(1, 3)]), as.vector(m[, c(2,
                4)]), type = "n", xlab = xlabs, ylab = ylabs,
                main = opz$main, sub = opz$sub, ylim = range(fit, rangeCI))
        if(rug) {segments(xvalues, rep(par()$usr[3],length(xvalues)), xvalues,
            rep(par()$usr[3],length(xvalues))+ abs(diff(par()$usr[3:4]))/40)}
            }
        if(conf.level>0){
          if(rev.sgn) vall<- -vall
          if(shade) polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])),
            col = "gray", border=NA) else matlines(vall, ciValues[,-1], type="l", lty=2, col=cols)
            }
        if (res) points(xvalues, fit, cex = cexs, pch = pchs, col = res.col)
        yhat <- x$family$linkinv(yhat)
        if (length(cols) == 1)
            cols <- rep(cols, max(id.group))
        if (length(lwds) == 1)
            lwds <- rep(lwds, max(id.group))
        if (length(ltys) == 1)
            ltys <- rep(ltys, max(id.group))
        for (i in 1:max(id.group)) {
            lines(xhat[id.group == i], yhat[id.group == i], col = cols[i],
                lwd = lwds[i], lty = ltys[i])
        }
    } else {
        r <- cbind(val, y.val)
        r1 <- cbind(val, y.val1)
        rr <- rbind(r, r1)
        fit <- c(y.val, y.val1)
        if (res) {
            ress <- if (inherits(x, what = "glm", which = FALSE))
                residuals(x, "working") * sqrt(x$weights)
            else resid(x)
            #fit <- broken.line(x, term, gap = show.gap, link = link, interc = TRUE) + ress + const
            fit <- predict.segmented(x, ifelse(rep(rev.sgn, length(xvalues)),-xvalues,xvalues), type=tipo) + ress + const
        }
        if (!add)
            plot(rr, type = "n", xlab = xlabs, ylab = ylabs,
                main = opz$main, sub = opz$sub, ylim = range(fit, rangeCI))
        if(rug) {segments(xvalues, rep(par()$usr[3],length(xvalues)), xvalues,
            rep(par()$usr[3],length(xvalues))+ abs(diff(par()$usr[3:4]))/40)}
        if(conf.level>0) {
          if(rev.sgn) vall<- -vall
          if(shade) polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])),
            col = "gray", border=NA) else matlines(vall, ciValues[,-1], type="l", lty=2, col=cols)
            }
        if (res)
            points(xvalues, fit, cex = cexs, pch = pchs, col = res.col)
            segments(m[, 1], m[, 2], m[, 3], m[, 4], col = cols, lwd = lwds, lty = ltys)
        }
    invisible(NULL)
}
