plot.stepmented <- function(x, term, add = FALSE, res = TRUE, conf.level=0, interc = TRUE, add.fx=FALSE,
                            psi.lines = TRUE, link = TRUE, const=NULL, res.col = grey(.15, alpha = .4), 
                            surf=FALSE, zero.cor=TRUE, heurs=TRUE, shade=FALSE, se.type=c("cdf","abs","none"), 
                            k=NULL, .vcov=NULL, leg="topleft", ...) {
  #=============================
  plot2step<-function(object, res, interc, psi.lines, arg, add){
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
    estpsi1=estpsi[nomiZ[1]==gsub("psi[1234567890]*[1234567890].","",nomiPsi)]
    estpsi2=estpsi[nomiZ[2]==gsub("psi[1234567890]*[1234567890].","",nomiPsi)]
    
    #browser()
    if(is.null(arg$xlim)) arg$xlim<-object$rangeZ[,nomiZ[1]]
    if(is.null(arg$ylim)) arg$ylim<-object$rangeZ[,nomiZ[2]]
    
    #browser()
    
    if(!add){
      plot(1, xlab=nomiZ[1], ylab=nomiZ[2], xaxt="n", yaxt="n", ylim=arg$ylim, xlim=arg$xlim, xaxs="i", yaxs="i")
      axis(2, at= round( x$rangeZ[, nomiZ[2]]+c(diff(x$rangeZ[, nomiZ[2]])/100,0), 3) , cex.axis=.7)
      axis(2, at= round(c(estpsi2), 3), cex.axis=.7, las=2)
      axis(1, at= round(c(estpsi1, object$rangeZ[, nomiZ[1]]), 3), cex.axis=.7)
    }

    xvalues<-sort(c(par()$usr[1:2], estpsi1))
    yvalues<-sort(c(par()$usr[3:4], estpsi2))
    
    #browser()
    fcol<-function(.x,R=range(.x))colorRampPalette(c(arg$col, "white"), alpha=TRUE)(200)[findInterval(.x,seq(min(R),max(R),l=200))]    

    #fit01<-fit01[nrow(fit01):1,]
    for(j in 1:(length(estpsi1)+1)){
      A=cbind(rep(xvalues[j],length(yvalues)-1), yvalues[-length(yvalues)])
      B=cbind(rep(xvalues[j+1],length(yvalues)-1), yvalues[-1])
      cc=if(arg$col==1) grey(fit01[,j], alpha=.75) else fcol(fit01[,j], range(fit01)) 
      if(add) cc<-adjustcolor(cc, alpha.f=.5)
      rect(A[,1],A[,2],B[,1],B[,2], col=cc, border=NA)
    }
    
    
    xx<-xvalues[-length(xvalues)]+diff(xvalues)/2
    yy<- yvalues[-length(yvalues)]+diff(yvalues)/2
    
    E<-expand.grid(yy,xx)
    
    #browser()
    if(!add){
    if(res) {
      r<-object$residuals+object$fitted.values
      r<- (r-min(r))/diff(range(r))
      r<-arg$cex/2+ (arg$cex*3-arg$cex/2)*r 
      points(object$Z, pch=arg$pch, col=adjustcolor(arg$res.col,.5), cex=r)
      #points(object$Z, pch=arg$pch, cex=arg$cex, col=adjustcolor(arg$res.col,.5))
    }
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
  pred.step.plot<-function(object, k=NULL, apprx=c("cdf","abs","none"), nomeZ, zero.cor=TRUE, V=NULL, ...){
    apprx=match.arg(apprx)
    X=model.matrix.stepmented(object, k=k, type=apprx)
    if(is.null(V)) V<-vcov.stepmented(object, zero.cor=zero.cor, type=apprx)

    interc=TRUE
    nomiZ<- object$nameUV$Z
    nomiV<- object$nameUV$V
    nomiU<- object$nameUV$U
    nomiPsi<- gsub("V","psi", nomiV)
    #id.noV<-setdiff(colnames(X), nomiPsi)
    
    #browser()
    
    #se ci sono piu' variabili segmented seleziona solo i termini di un unica variabile
    #if(missing(nomeZ)) nomeZ <-nomiZ[1]
    nomeZ <- if(missing(nomeZ)) nomiZ[1] else unique(nomeZ)
    nomeU.ok   <- grep(paste(".", nomeZ,sep=""), nomiU, value=TRUE)
    nomePsi.ok <- grep(paste(".", nomeZ,sep=""), nomiPsi, value=TRUE)
    
    N=20
    vv<-matrix(, N, length(nomePsi.ok))
    for(i in 1:length(nomePsi.ok)){
      #psi <- object$psi[nomePsi.ok[i],"Est."]
      psi <- object$psi.rounded[1, nomePsi.ok[i]]
      #se  <- object$psi[nomePsi.ok[i],"St.Err"]
      se  <- sqrt(diag(V))[nomePsi.ok[i]]
      if(is.na(se)) se<-1e-5
      vv[,i] <-  qnorm(seq(.0001, .9999,l=N), psi, se) #mettere 2*se???
    }
    vv<-as.vector(vv)  
    mi=object$rangeZ[1,nomeZ]
    ma=object$rangeZ[2, nomeZ]
    vv<-sort(c(object$psi.rounded[1, nomePsi.ok], seq(mi,ma,l=30), vv))
    vv<- vv[vv>=mi & vv<= ma]
    
    #browser()
    
    if(apprx!="none"){
      Xok<-matrix(, length(vv), 2*length(nomePsi.ok))
      colnames(Xok)<- c(nomeU.ok,nomePsi.ok)
      for(i in 1:length(nomePsi.ok)){
        Xok[, nomeU.ok[i]]   <- spline(object$Z[,nomeZ], X[, nomeU.ok[i]], xout=vv)$y
        Xok[, nomePsi.ok[i]] <- spline(object$Z[,nomeZ], X[, nomePsi.ok[i]], xout=vv)$y
      }
    } else {
      Xok<-matrix(, length(vv), length(nomeU.ok))
      colnames(Xok)<- nomeU.ok
      for(i in 1:length(nomeU.ok)){
        Xok[, nomeU.ok[i]]   <- spline(object$Z[,nomeZ], X[, nomeU.ok[i]], xout=vv)$y
        #Xok[, nomePsi.ok[i]] <- spline(object$Z[,nomeZ], X[, nomePsi.ok[i]], xout=vv)$y
      }
    }
    
    id.ok = grep(paste(".", nomeZ,sep=""), colnames(X))
    
    psii = object$psi.rounded[1, nomePsi.ok]
    #psii = object$psi[nomePsi.ok,"Est."]
    M=sapply(psii, function(.x) 1*(vv>.x))
    cof<-object$coefficients[nomeU.ok]
    id.interc = FALSE
    if(interc && "(Intercept)" %in% colnames(X)) {
      id.interc=TRUE
      id.ok<-c(1, id.ok)
      Xok<-cbind(1, Xok)
      M = cbind(1, M)
      cof<-c(object$coefficients[1], cof)
    }
    
    #browser()
    
    se=rowSums((Xok%*%V[id.ok,id.ok])*Xok)
    if(any(se<=0)) {
      warning("correcting non-positive st.err 1", call.=FALSE)
      se[se<=0]<-median(se[se>0])
    }
    se.smooth=sqrt(se)
    fit.smooth <- drop(Xok[,-match(nomePsi.ok, colnames(Xok), 0)]%*%cof)
    
    #=======================================
    Xok[, nomeU.ok] <- if(id.interc) M[,-1] else M
    se=rowSums((Xok%*%V[id.ok,id.ok])*Xok)
    if(any(se<=0)) {
      warning("correcting non-positive st.err 2", call.=FALSE)
      se[se<=0]<-median(se[se>0])
    }
    se.nonsmooth=sqrt(se)
    fit.nonsmooth<-drop(M %*% cof)
#browser()

    g <- cut(vv, breaks = c(mi, psii, ma), labels =FALSE, include.lowest = TRUE)
    r<-list(values=vv, g=g, fit.nonsmooth=fit.nonsmooth, se.nonsmooth=se.nonsmooth,
            fit.smooth=fit.smooth, se.smooth=se.smooth)
    r
  }
  #=============================
  
  
  arg <- list(...)
  se.type <- match.arg(se.type)
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
    plot2step(x, res=res, interc=interc, psi.lines=psi.lines, arg=arg, add=add)
    return(invisible(NULL))
  }
  
  #browser()
  
  if (missing(term)) {
    if (length(unique(x$nameUV$Z)) > 1) {
      stop("please, specify `term'")
    } else {
      term <- x$nameUV$Z
    }
  } else {
    if(is.numeric(term)) term <- x$nameUV$Z[term]
    #if(length(term)>=2) stop(" 'term' should be a scalar")
    #dterm <- deparse(substitute(term))
    #if (dterm %in% x$nameUV$Z) term <- dterm
    if (!isTRUE(all(term %in% x$nameUV$Z))) stop(paste("Unknown term. It should be numeric or one of: ", paste(" '", x$nameUV$Z, "' ", sep="", collapse="")))
  }
  if(add.fx && !term%in%colnames(x$f.x)) stop("no additional effect for the selected term")
  
  if(is.null(.vcov)) .vcov <- vcov.stepmented(x, zero.cor=zero.cor, type=se.type) 
  
  
  if(length(term)>1){
    opz<-list(...)
    cols<- if(!is.null(opz$col)) opz$col else 1:length(term)+1
    cols <- rep(cols, l=length(term)) 
    res.cols<- rep(res.col, l=length(term))
    lwds<- if("lwd"%in% names(opz)) opz$lwd else 2
    lwds<- rep(lwds, l=length(term))
    ltys<- if("lty"%in% names(opz)) opz$lty else 1
    ltys<- rep(ltys, l=length(term))
    cexs<- if("cex"%in% names(opz)) opz$cex else .75
    cexs<- rep(cexs, l=length(term))
    pchs<- if("pch"%in% names(opz)) opz$pch else 19
    pchs<- rep(pchs, l=length(term))
    if(!is.null(opz$ylim)) {
      Ylim <- opz$ylim 
    } else {
      if(inherits(x, "glm")){
        if(link){
          Ylim <- if(!res) range(x$linear.predictors) else range(x$linear.predictors+x$residuals)
        } else {
          Ylim <- if(!res) range(x$fitted.values) else range(x$fitted.values+ residuals(x, "response"))
        }
      } else {
        Ylim <- if(!res) range(x$fitted.values) else range(x$fitted.values+x$residuals)
      }
      
    }
    Ylab <- if(!is.null(opz$ylab)) opz$ylab else paste(formula(x))[2]
    idTerm <- if(is.numeric(term)) term else match(term, x$nameUV$Z)
    nomeX <- intersect(strsplit(x$nameUV$Z,":")[[idTerm[1]]], unlist(strsplit(x$nameUV$Z,":")[idTerm[-1]]))
    Xlab <- if(!is.null(opz$xlab)) opz$xlab else nomeX
    Xlim<- if(!is.null(opz$xlim)) opz$xlim else range(x$rangeZ[,term]) #range(x$model[,nomeX])
    int.all<-rep(NA, length(term))
    
    plot.stepmented(x, term[1], add = add, res = res, conf.level = conf.level, .vcov=.vcov,
                   interc=interc, link = link, res.col = res.cols[1], const = const, shade=shade,
                    surf=FALSE, zero.cor = zero.cor, add.fx=add.fx, heurs=heurs,  se.type=se.type, k=k, psi.lines=psi.lines,
                   #rug=FALSE, dens.rug=FALSE, dens.col = grey(0.8), smoos=NULL, hide.zeros=TRUE, 
                   #transf=I, isV=FALSE, is=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, .coef=NULL, prev.trend=FALSE, 
                   col=cols[1], ylim=Ylim, xlim=Xlim, ylab=Ylab,xlab=Xlab,
                   lty=ltys[1],pch=pchs[1],lwd=lwds[1],cex=cexs[1])
    Term<- if(is.numeric(term[1])) x$nameUV$Z[term[1]] else term[1]
    int.all[1]<-interc.gr<- strsplit(Term, ":")[[1]][2]
    #points.segmented(x, term[1], col=cols[1], const=estcoef[interc.gr], v=psi.lines, pch=20, link=link)
    for(j in 2:length(term)){
      #browser()
      plot.stepmented(x, term[j], add = TRUE, res = res, conf.level = conf.level, .vcov=.vcov,
                      interc=interc, link = link, res.col = res.cols[j], const = const, shade=shade,
                      surf=FALSE, zero.cor = zero.cor, add.fx=add.fx, heurs=heurs,  se.type=se.type, k=k, psi.lines=psi.lines,
                      #rug=FALSE, dens.rug=FALSE, dens.col = grey(0.8), smoos=NULL, hide.zeros=TRUE, 
                      #transf=I, isV=FALSE, is=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, .coef=NULL, prev.trend=FALSE, 
                      col=cols[j], ylim=Ylim, xlim=Xlim, ylab=Ylab,xlab=Xlab,
                      lty=ltys[j], pch=pchs[j],lwd=lwds[j], cex=cexs[j]) #hide.zeros smoos=NULL
      
      # plot.stepmented(x, term[j], add = TRUE, res = res, conf.level = conf.level, 
      #                interc=interc, link = link, res.col = res.cols[j], rev.sgn = rev.sgn, const = const, 
      #                shade=shade, rug=FALSE, dens.rug=FALSE, dens.col = grey(0.8),
      #                transf=I, isV=FALSE, is=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, .coef=NULL, prev.trend=FALSE, 
      #                smoos=NULL, hide.zeros=TRUE,col=cols[j],
      #                lty=ltys[j],pch=pchs[j],lwd=lwds[j],cex=cexs[j])
      Term<- if(is.numeric(term[j])) x$nameUV$Z[term[j]] else term[j]
      int.all[j]<-interc.gr<- strsplit(Term, ":")[[1]][2]
      #points.segmented(x, term[j], col=cols[j], const = estcoef[interc.gr], v=psi.lines, pch=20, link=link)
    }
    if(!is.na(leg)) {
      legend(leg, int.all, col=cols, lty=1, lwd=1.5, bty="n")
    }
  } else { #se length(term)==1
    estcoef <- x$coefficients
    if(is.null(const)){
      interc.gr<- strsplit(term, ":")[[1]][2]
      const<- estcoef[interc.gr]
      if(is.na(const)) const<-0
    }
    #browser()
  
    idU <- x$nameUV$U[endsWith(x$nameUV$U, paste(".", term, sep = ""))]
    if (interc && "(Intercept)" %in% names(x$coefficients)) {
      idU <- c("(Intercept)", idU)
      est.means <- cumsum(x$coefficients[idU])
    } else {
      est.means <- c(0,cumsum(x$coefficients[idU]))
    }
    
    nomiPsi <- sub("V", "psi", x$nameUV$V)
    idPsi <- nomiPsi[endsWith(nomiPsi, paste(".", term, sep = ""))]
    #psi <- sort(x$coefficients[idPsi])
    psi <- sort(x$psi.rounded[1,idPsi])
    rangeZ <- x$rangeZ[, term]
    Z <- drop(x$Z[, term, drop = TRUE])
    m <- min(rangeZ)
    M <- max(rangeZ)
    #browser()
    
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
  
    y<-y+const
    Y<-Y+const
    
    
    if (inherits(x, what = "glm", which = FALSE) && !link){
      y<- x$family$linkinv(y)
      Y<- x$family$linkinv(Y)
      if(add.fx && !is.null(x$f.x)) f.x100<- x$family$linkinv(f.x100)
    }
    
    if(conf.level>0){
      if(add.fx && !is.null(x$f.x)) stop(" conf.level>0 is not allowed with additional terms")
      r<-pred.step.plot(x, k=k, apprx=se.type, nomeZ=term, zero.cor=zero.cor, V=.vcov)
      #browser()
      vv<-r$values
      ff<-r$fit.nonsmooth
      se<-r$se.nonsmooth
      z<-if(inherits(x, "glm")) abs(qnorm((1-conf.level)/2)) else abs(qt((1-conf.level)/2, df= x$df.residual))
      inf= ff-z*se
      sup= ff+z*se
      if(heurs && se.type!="none") {
        inf<-pmin(inf, r$fit.smooth-z*r$se.smooth )
        sup<-pmax(sup, r$fit.smooth+z*r$se.smooth )
      }
      #ff=r$fit.smooth
      #inf=r$fit.smooth-z*r$se.smooth
      #sup=r$fit.smooth+z*r$se.smooth
      M=cbind(vv, ff, inf, sup)
      if(is.null(arg$ylim)) arg$ylim<-range(c(inf, sup))
    } 
    
    return.ci=FALSE #era per fare alcune verifiche..
    if(return.ci){
      return(M)
    } else {
      ###############################################################  
      if(res) y<- if(link) y+x$residuals else y + resid(x, "response")
      if(add) {
        if(res) points(Z, y, pch = arg$pch, col =res.col, cex = arg$cex)
        } else {
          if (res) {
            plot(Z, y, ylab = arg$ylab, xlab = arg$xlab, pch = arg$pch, col =res.col,
              cex = arg$cex, xlim = arg$xlim, ylim = arg$ylim,
              main=arg$main, sub=arg$sub, cex.axis=arg$cex.axis, cex.lab=arg$cex.lab)
            } else {
              plot(Z, y, ylab = arg$ylab, xlab = arg$xlab, type = "n",
                xlim = arg$xlim, ylim = arg$ylim,
                main=arg$main, sub=arg$sub, cex.axis=arg$cex.axis, cex.lab=arg$cex.lab)
            }
        }
      coll <- rep(arg$col, length(est.means))
      ltyy <- rep(arg$lty, length(est.means))
      lwdd <- rep(arg$lwd, length(est.means))
      
      
      
      if(add.fx && !is.null(x$f.x)){
        #limy<-tapply(f.x100, g, function(.x).x[1])
        limy<-tapply(f.x100, g, function(.x)c(.x[1], .x[length(.x)]))
        limy<- unlist(limy)
        limy<-c(limy[1], rep(limy[2:(length(limy)-1)], each=2), limy[length(limy)])
        limy<- matrix(limy, ncol=2, byrow=TRUE)
        limy<- apply(limy[-c(1, nrow(limy)),,drop=FALSE],1,max)
        #limy<- do.call(rbind, limy)
        for(i in 1:length(est.means)) lines(Z100a[g==i], f.x100[g==i],
                                          col = coll[i], lwd = lwdd[i], lty = ltyy[i])
                                    #col = arg$col, lwd = arg$lwd, lty = arg$lty)
        } else {
          limy<- apply(matrix(Y[-c(1,length(Y))], ncol = 2, byrow = TRUE), 1, max)
          segments(X[id1], Y[id1], X[id2], Y[id2],
               col = coll, lwd = lwdd, lty = ltyy)
               #col = arg$col, lwd = arg$lwd, lty = arg$lty)
        }
        # ########################################################
    }
    
     
    #browser()
    
    if(conf.level>0){
      #browser()
      if(shade) {
        polygon(c(vv, rev(vv)), c(M[,3],rev(M[,4])),
                col = adjustcolor(coll, .3), border=NA) 
      } else {
        for(i in 1:length(est.means)){
          matlines(vv[r$g==i], M[r$g==i, c(2,3,4)], lty=c(1,2,2), col=coll[i])
        }
      }
      #matlines(vv, cbind(ff, ff-z*se, ff+z*se), lty=c(1,2,2), col=coll)  
    }
    
    # plot(X,Y, type='s)
    #browser()
    if(psi.lines) {
      #limy<-apply(matrix(Y[-c(1,length(Y))], ncol = 2, byrow = TRUE), 1, max)
      segments(x0 = psi, y0 = par()$usr[3], x1 = psi, y1 = limy, lty = 3, col = 1) #arg$col)
      points(psi, rep(par()$usr[3], length(psi)), pch = 19, col = arg$col)
    }
    invisible(NULL)
  }
}
