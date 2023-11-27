plot.stepmented <- function(x, term, add = FALSE, res = TRUE, interc = TRUE, add.fx=FALSE,
                                  psi.lines = TRUE, link = FALSE, surf=FALSE, ...) {
  #=============================
  plot2step<-function(object, res, interc, psi.lines, arg){
    nomiZ<- object$nameUV$Z
    if(length(nomiZ)!=2) stop("surf=TRUE is allowed only with 2 stepmented covariates")
    nomiV<- object$nameUV$V
    nomiPsi<-gsub("V","psi",nomiV)
    nomiU<- object$nameUV$U
    
    m1<-slope(object)[[nomiZ[1]]][,"Est."]-coef(object)["(Intercept)"]
    m2<-slope(object)[[nomiZ[2]]][,"Est."]
    if(!interc){
      m1<-m1 
      m2<-m2 -coef(object)["(Intercept)"]
    }
    fit <-outer(m2,m1,"+")
    fit01<- 1- (fit-min(fit))/diff(range(fit))
    
    estpsi=object$psi[,"Est."]
    estpsi1=estpsi[nomiZ[1]==gsub("psi[1234567890].","",nomiPsi)]
    estpsi2=estpsi[nomiZ[2]==gsub("psi[1234567890].","",nomiPsi)]
    
    #browser()
    if(is.null(arg$xlim)) arg$xlim<-object$rangeZ[,nomiZ[1]]
    if(is.null(arg$ylim)) arg$ylim<-object$rangeZ[,nomiZ[2]]
    plot(1, xlab=nomiZ[1], ylab=nomiZ[2], xaxt="n", yaxt="n", ylim=arg$ylim, xlim=arg$xlim, xaxs="i", yaxs="i")
    
    #browser()
    axis(2, at= round( x$rangeZ[, nomiZ[2]]+c(diff(x$rangeZ[, nomiZ[2]])/100,0), 3) , cex.axis=.7)
    axis(2, at= round(c(estpsi2), 3), cex.axis=.7, las=2)
    axis(1, at= round(c(estpsi1, object$rangeZ[, nomiZ[1]]), 3), cex.axis=.7)
    
    xvalues<-sort(c(par()$usr[1:2], estpsi1))
    yvalues<-sort(c(par()$usr[3:4], estpsi2))
    
    #browser()
    fcol<-function(.x,R=range(.x))colorRampPalette(c(arg$col, "white"), alpha=TRUE)(200)[findInterval(.x,seq(min(R),max(R),l=200))]    

    #fit01<-fit01[nrow(fit01):1,]
    for(j in 1:(length(estpsi1)+1)){
      A=cbind(rep(xvalues[j],length(yvalues)-1), yvalues[-length(yvalues)])
      B=cbind(rep(xvalues[j+1],length(yvalues)-1), yvalues[-1])
      cc=if(arg$col==1) grey(fit01[,j], alpha=.75) else fcol(fit01[,j], range(fit01)) 
      rect(A[,1],A[,2],B[,1],B[,2], col=cc, border=NA)
    }
    
    
    xx<-xvalues[-length(xvalues)]+diff(xvalues)/2
    yy<- yvalues[-length(yvalues)]+diff(yvalues)/2
    
    E<-expand.grid(yy,xx)
    
    #browser()
    
    if(res) {
      r<-object$residuals+object$fitted.values
      r<- (r-min(r))/diff(range(r))
      r<-arg$cex/2+ (arg$cex*3-arg$cex/2)*r 
      points(object$Z, pch=arg$pch, col=adjustcolor(arg$res.col,.5), cex=r)
      #points(object$Z, pch=arg$pch, cex=arg$cex, col=adjustcolor(arg$res.col,.5))
    }
    
    if(psi.lines) abline(v=estpsi1, h=estpsi2, lty=arg$lty)
    
    #browser()
    
    #text(E[,2],E[,1], round(as.vector(fit[nrow(fit):1,]),3))
    text(E[,2],E[,1], round(as.vector(fit),3))
    
    box()
    return(invisible(NULL))
    #X<-matrix(,, nrow(o$psi))
    #X[,i]<- -(model.matrix(o)[,nomiU[i]]-.5)/model.matrix(o)[,nomiPsi[i]]
    #nomiPsi nomiU
  }  
  
  #=============================
  arg <- list(...)
  if (is.null(arg$col))    arg$col = 2#grey(0.4)
  if (is.null(arg$lwd))    arg$lwd = 2.5
  if (is.null(arg$lty))    arg$lty = 1
  
  if (is.null(arg$res.col))arg$res.col= grey(0.7, alpha=.8)
  if (is.null(arg$cex))    arg$cex = 1.2
  if (is.null(arg$pch))    arg$pch = 20

  if (is.null(arg$xlim))   arg$xlim = NULL
  if (is.null(arg$ylim))   arg$ylim = NULL
  if (is.null(arg$main))   arg$main = NULL
  if (is.null(arg$sub))   arg$sub = NULL
  if (is.null(arg$cex.axis))   arg$cex.axis = 1
  if (is.null(arg$cex.lab))   arg$cex.lab = 1

  if(surf){
    if(length(x$nameUV$Z)!=2) stop(" 'surf=TRUE' works only with 2 stepmented terms")
    plot2step(x, res=res, interc=interc, psi.lines=psi.lines, arg=arg)
    return(invisible(NULL))
  }
  
  if (missing(term)) {
    if (length(x$nameUV$Z) > 1) {
      stop("please, specify `term'")
    } else {
      term <- x$nameUV$Z
    }
  } else {
    dterm <- deparse(substitute(term))
    if (dterm %in% x$nameUV$Z)
      term <- dterm
    if (!isTRUE(term %in% x$nameUV$Z))
      stop("invalid `term'")
  }
  if(add.fx && !term%in%colnames(x$f.x)) stop("no additional effect for the selected term")
  idU <- x$nameUV$U[endsWith(x$nameUV$U, paste(".", term, sep = ""))]
  if (interc && "(Intercept)" %in% names(x$coefficients)) idU <- c("(Intercept)", idU)
  est.means <- cumsum(x$coefficients[idU])
  nomiPsi <- gsub("V", "psi", x$nameUV$V)
  idPsi <- nomiPsi[endsWith(nomiPsi, paste(".", term, sep = ""))]
  psi <- sort(x$coefficients[idPsi])
  rangeZ <- x$rangeZ[, term]
  Z <- drop(x$Z[, term, drop = TRUE])
  m <- min(rangeZ)
  M <- max(rangeZ)
  fit0 <- as.numeric(as.character(cut(Z, breaks = c(m, psi, M), labels = est.means, include.lowest = TRUE)))
  #y <- x$residuals + fit0  #oppure: x$model[[1]]
  
  #browser()
  
  y<- fit0
  
  Y <- rep(est.means, each = 2)
  X <- c(m, rep(psi, each = 2), M)
  id1 <- seq(1, length(X), by = 2)
  id2 <- seq(2, length(X), by = 2)
  
  if (is.null(arg$xlab)) arg$xlab = colnames(x$Z[, term, drop = FALSE])
  if (is.null(arg$ylab)) arg$ylab = names(x$model)[1]
  
  #browser()
  #se c'e' un effetto della x aggiuntivo
  if(add.fx && !is.null(x$f.x)){
    Z100<- seq(min(Z), max(Z), l=nrow(x$f.x))
    #Z100a<- sort(c(psi*c(.999,1,1.0001),seq(min(Z), max(Z), l=nrow(x$f.x)))) #-3*length(psi))))
    Z100a<- sort(c(as.vector(sapply(psi, function(.x) .x*c(.999,1,1.0001))), seq(min(Z), max(Z), l=nrow(x$f.x))))
    g <- cut(Z100a, breaks = c(m, psi, M), labels =FALSE, include.lowest = TRUE)
    f.x100 <- x$f.x[, term, drop=TRUE]
    fSpline<-splinefun(Z100, f.x100)
    f.x.n <-fSpline(Z)
    y<- y + f.x.n
    f.x100 <- fSpline(Z100a)+ rep(est.means, as.numeric(table(g))) #cbind(Z100a, g, rep(est.means, as.numeric(table(g))))
    #fit0<- fit0 + f.x.n 
    #y<- y + x$f.x[,term, drop=TRUE]
    #fit0<- fit0 + x$f.x[,term, drop=TRUE]
  }

  #browser()
  if (inherits(x, what = "glm", which = FALSE) && !link){
    y<- x$family$linkinv(y)
    Y<- x$family$linkinv(Y)
    if(add.fx && !is.null(x$f.x)) f.x100<- x$family$linkinv(f.x100)
  }
  
  if(res) {
    y<- if(link) y+x$residuals else y + resid(x, "response")
  }
  
  if(!add){
    if (res) {
      plot(Z, y, ylab = arg$ylab, xlab = arg$xlab, pch = arg$pch, col =arg$res.col,
        cex = arg$cex, xlim = arg$xlim, ylim = arg$ylim,
        main=arg$main, sub=arg$sub, cex.axis=arg$cex.axis, cex.lab=arg$cex.lab)
      } else { 
      plot(Z, y, ylab = arg$ylab, xlab = arg$xlab, type = "n", 
        xlim = arg$xlim, ylim = arg$ylim, 
        main=arg$main, sub=arg$sub, cex.axis=arg$cex.axis, cex.lab=arg$cex.lab)
      }
  }
   coll <- rep(arg$col, length(est.means))
   ltyy <-rep(arg$lty, length(est.means))
   lwdd <- rep(arg$lwd, length(est.means))
  #browser()
   if(add.fx && !is.null(x$f.x)){
    limy<-tapply(f.x100, g, function(.x).x[1])
    for(i in 1:length(est.means)) lines(Z100a[g==i], f.x100[g==i],
                                        col = coll[i], lwd = lwdd[i], lty = ltyy[i])
                                  #col = arg$col, lwd = arg$lwd, lty = arg$lty)
    
      } else {
    limy<- apply(matrix(Y[-c(1,length(Y))], ncol = 2, byrow = TRUE), 1, max)
    segments(X[id1], Y[id1], X[id2], Y[id2], 
             col = coll, lwd = lwdd, lty = ltyy)
             #col = arg$col, lwd = arg$lwd, lty = arg$lty)
  }
  # plot(X,Y, type='s)
  if(psi.lines) {
    #limy<-apply(matrix(Y[-c(1,length(Y))], ncol = 2, byrow = TRUE), 1, max)
    segments(x0 = psi, y0 = par()$usr[3], x1 = psi, y1 = limy, lty = 3, col = 1) #arg$col)
    points(psi, rep(par()$usr[3], length(psi)), pch = 19, col = arg$col)
  }
  invisible(NULL)
}
