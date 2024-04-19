plot.segmented<-function (x, term, add = FALSE, res = FALSE, conf.level = 0, 
                          interc=TRUE, link = TRUE, res.col = grey(.15, alpha = .4), rev.sgn = FALSE, const = NULL, 
                          shade=FALSE, rug=!add, dens.rug=FALSE, dens.col = grey(0.8),
                          transf=I, isV=FALSE, is=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, .coef=NULL, prev.trend=FALSE, 
                          smoos=NULL, hide.zeros=FALSE, leg="topleft", psi.lines=FALSE, ...){
  #put leg=NA if you do not want the legend..
  #funzione plot.segmented che consente di disegnare anche i pointwise CI
  f.U<-function(nomiU, term=NULL){
    #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
    #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
    k<-length(nomiU)
    nomiUsenzaU<-strsplit(nomiU, "\\.")
    nomiU.ok<-vector(length=k)
    for(i in 1:k){
      nomi.i<-nomiUsenzaU[[i]][-1]
      if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
      nomiU.ok[i]<-nomi.i
    }
    if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
    return(nomiU.ok)
  }
  #-------------- 
  enl.range<-function(..., enlarge=TRUE){
    #modifica il min dei valori in ...
    r<-range(...)
    if(enlarge) r[1]<-if(sign(r[1])>0) r[1]*.9 else r[1]*1.1
    r
  }
  #--------------
  #se l'oggetto e' segmented.Arima il nome dell'eventuale interc va sostituito..
  #if((all(class(x)==c("segmented", "Arima")))) names(x$coef)<-gsub("intercept", "(Intercept)", names(coef(x)))
  if(all(c("segmented", "Arima") %in% class(x))) names(x$coef)<-gsub("intercept", "(Intercept)", names(x$coef))
  
  covv <- if(is.null(.vcov)) vcov(x, is=is, var.diff=var.diff) else .vcov 

  if(!is.null(.coef)) {
    estcoef<- .coef
    } else { 
      estcoef <- coef(x)
      if(is.null(estcoef)) estcoef <- x$coef
      if(is.null(estcoef)) stop("No coeffs in the fit? Please use '.coef'")
    }
  
  if(length(estcoef)==0) stop("No coefficient in the object fit?")
  
  #browser()
  if(!all(dim(covv)==c(length(estcoef), length(estcoef)))) stop("dimension of cov matrix and estimated coeffs do not match", call. = FALSE)
  
  #--------------
  linkinv <- !link
  if (inherits(x, what = "glm", which = FALSE) && linkinv && !is.null(x$offset) && res) stop("residuals with offset on the response scale?")
  if(conf.level< 0 || conf.level>.9999) stop("meaningless 'conf.level'")
  if ((inherits(x, what = "glm", which = FALSE) && linkinv) || res) {
    if(!(identical(transf, I) || identical(transf, "I"))) {transf<-I; warning("'transf' set to I with 'res=TRUE' and/or 'link=FALSE'.")}
  }
  if(missing(term)) {
    if (length(x$nameUV$Z) > 1) {
      stop("please, specify `term'")
    } else {
      term <- x$nameUV$Z
    }
  } else {
    #browser()
    if(is.numeric(term)) term <- x$nameUV$Z[term]
    #if(!is.character(term)) stop("please specify correctly 'term' ")
    #term<- deparse(substitute(term))
    #if(dterm %in% x$nameUV$Z) term<-dterm
    if (!isTRUE(all(term %in% x$nameUV$Z))) stop(paste("Unknown term. It should be numeric or one of: ", paste(" '", x$nameUV$Z, "' ", sep="", collapse="")))
  }
  
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
    Ylim <- if(!is.null(opz$ylim)) opz$ylim else range(x$fitted.values+x$residuals)
    Ylab <- if(!is.null(opz$ylab)) opz$ylab else paste(formula(x))[2]
    idTerm <- if(is.numeric(term)) term else match(term, x$nameUV$Z)
    nomeX <- intersect(strsplit(x$nameUV$Z,":")[[idTerm[1]]], unlist(strsplit(x$nameUV$Z,":")[idTerm[-1]]))
    Xlab <- if(!is.null(opz$xlab)) opz$xlab else nomeX
    Xlim<- if(!is.null(opz$xlim)) opz$xlim else range(x$model[,nomeX])
    int.all<-rep(NA, length(term))

    plot.segmented(x, term[1], add = add, res = res, conf.level = conf.level, 
                   interc=interc, link = link, res.col = res.cols[1], rev.sgn = rev.sgn, const = const, 
                   shade=shade, rug=FALSE, dens.rug=FALSE, dens.col = grey(0.8),
                   transf=I, isV=FALSE, is=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, .coef=NULL, prev.trend=FALSE, 
                   smoos=NULL, hide.zeros=TRUE, col=cols[1], ylim=Ylim, xlim=Xlim, ylab=Ylab,xlab=Xlab,
                   lty=ltys[1],pch=pchs[1],lwd=lwds[1],cex=cexs[1])
    Term<- if(is.numeric(term[1])) x$nameUV$Z[term[1]] else term[1]
    int.all[1]<-interc.gr<- strsplit(Term, ":")[[1]][2]
    points.segmented(x, term[1], col=cols[1], const=estcoef[interc.gr], v=psi.lines, pch=20)
    for(j in 2:length(term)){
      plot.segmented(x, term[j], add = TRUE, res = res, conf.level = conf.level, 
                     interc=interc, link = link, res.col = res.cols[j], rev.sgn = rev.sgn, const = const, 
                     shade=shade, rug=FALSE, dens.rug=FALSE, dens.col = grey(0.8),
                     transf=I, isV=FALSE, is=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, .coef=NULL, prev.trend=FALSE, 
                     smoos=NULL, hide.zeros=TRUE,col=cols[j],
                     lty=ltys[j],pch=pchs[j],lwd=lwds[j],cex=cexs[j])
      Term<- if(is.numeric(term[j])) x$nameUV$Z[term[j]] else term[j]
      int.all[j]<-interc.gr<- strsplit(Term, ":")[[1]][2]
      points.segmented(x, term[j], col=cols[j], const = estcoef[interc.gr], v=psi.lines, pch=20)
    }
    if(!is.na(leg)) {
      legend(leg, int.all, col=cols, lty=1, lwd=1.5, bty="n")
    }
  } else {
    if(is.null(const)){
      interc.gr<- strsplit(term, ":")[[1]][2]
      const<- estcoef[interc.gr]
      if(is.na(const)) const<-0
    }
    if(!is.numeric(const)) stop(" 'const' should be NULL (default) or numeric")
    opz <- list(...)
    col.shade<-if(!is.null(opz$col.shade)) adjustcolor(opz$col.shade, .15) else adjustcolor("grey", .4)
    cols<- if("col"%in% names(opz)) opz$col else 2
    lwds<- if("lwd"%in% names(opz)) opz$lwd else 2
    ltys<- if("lty"%in% names(opz)) opz$lty else 1
    cexs<- if("cex"%in% names(opz)) opz$cex else .75
    pchs<- if("pch"%in% names(opz)) opz$pch else 19
    ylabs<- if("ylab"%in% names(opz)) opz$ylab else paste("Effect  of ", term, sep = " ")
    xlabs<- if("xlab"%in% names(opz)) opz$xlab else term
  
    a <- intercept(x, term, digits=20, .vcov=covv, .coef=estcoef)[[1]][, "Est."]
    #Poiche' intercept() restituisce quantita' che includono sempre l'intercetta del modello, questa va eliminata se interc=FALSE
  
    idInterc<-grep("ntercept",names(estcoef))
    if(!interc && length(idInterc)==1) a<- a-estcoef[idInterc]
    b <- slope(x, term, digits=20, .coef=estcoef, .vcov=covv)[[1]][, "Est."]
  
    id <- f.U(rownames(x$psi), term)
    est.psi <- x$indexU[[term]]
    val <- sort(c(est.psi, x$rangeZ[, term]))
    #vettorializza i cols, lwds, ltys
    cols<-rep(cols, l=length(est.psi)+1)
    lwds<-rep(lwds, l=length(est.psi)+1)
    ltys<-rep(ltys, l=length(est.psi)+1)
    #---------aggiunta per gli IC
    rangeCI<-NULL
    vall<-sort(c(seq(min(val), max(val), l=100), est.psi, est.psi+1e-5))
    #ciValues<-predict.segmented(x, newdata=vall, se.fit=TRUE, type=tipo, level=conf.level)
    vall.list<-list(vall)
    names(vall.list)<-term
  
    if(conf.level>0) {
      k.alpha<- if(all(c("segmented","lm") %in% class(x))) abs(qt((1-conf.level)/2, x$df.residual)) else abs(qnorm((1-conf.level)/2))
      ciValues<-broken.line(x, vall.list, link=link, interc=interc, se.fit=TRUE, isV=isV, is=is, var.diff=var.diff, 
                          p.df=p.df, .vcov=covv, .coef=estcoef) #se gli passi covv, gli argomenti is e var.diff NON servono perche li ignora..
      ciValues<-cbind(ciValues$fit, ciValues$fit- k.alpha*ciValues$se.fit, ciValues$fit + k.alpha*ciValues$se.fit) + const
      #---> transf...
      ciValues<-apply(ciValues, 2, transf)
      rangeCI<-range(ciValues)
      #ciValues  e' una matrice di length(val)x3. Le 3 colonne: stime, inf, sup
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
  
    #xvalues <-  if(all(c("segmented", "Arima") %in% class(x))) x$Z[,1] else  model.matrix(x)[,term] #x$model[, term]
  
    #browser()
    if(inherits(x,"Arima")){
      xvalues <-x$Z[,1]
    } else {
      M <- model.matrix.segmented(x)
      #il 18/4/24 mi sono accorto che con ogg ottenuti da segmented.* con leftmost pendenza nulla non funzionava
      #perche' model.matrix.segmented non restituiva la variabile (non inserita nel modello (g)lm di partenza..)
      if(!term %in% colnames(M) && term%in%names(x$model)) M<-cbind(M, x$model[,term,drop=FALSE] )
      if(term %in% colnames(M)) {
        xvalues <- M[,term]
        } else {
          id.segTerm<-which(sapply(names(x$nameUV$formulaSeg), function(.x) startsWith(term,.x)))
          xvalues <- model.matrix(x$nameUV$formulaSeg[[id.segTerm]], data=x$model)[,term]
        }
    }

    if(rev.sgn) {
      val <- -val
      xvalues <- -xvalues
    }
    m <- cbind(val[s], y.val1[s], val[s + 1], y.val[s + 1])
  #values where to compute predictions (useful only if res=TRUE)
  #browser()
    if(res){
      new.d<-data.frame(ifelse(rep(rev.sgn, length(xvalues)),-xvalues, xvalues))
      names(new.d)<-term
      fit0 <- broken.line(x, new.d, link = link, interc=interc, se.fit=FALSE, .vcov=covv, .coef=estcoef)$fit
    }

  #-------------------------------------------------------------------------------

  if (inherits(x, what = "glm", which = FALSE) && linkinv) { #se GLM con link=FALSE (ovvero linkinv=TRUE)
    fit <- if (res)
      #predict.segmented(x, ifelse(rep(rev.sgn, length(xvalues)),-xvalues,xvalues), type=tipo) + resid(x, "response") + const
      #broken.line(x, term, gap = show.gap, link = link) + resid(x, "response") + const
      fit0 + resid(x, "response") + const        
    else x$family$linkinv(c(y.val, y.val1))
    xout <- sort(c(seq(val[1], val[length(val)], l = 50), val[-c(1, length(val))], 
                   pmax(val[-c(1, length(val))]*1.0001, val[-c(1, length(val))]*.9999)))
    l <- suppressWarnings(approx(as.vector(m[, c(1, 3)]), as.vector(m[, c(2, 4)]), xout = xout))
    val[length(val)]<- if(rev.sgn) min(l$x) else max(l$x) #aggiunto 11/09/17.. if else il 9/3/21
    
    id.group <- cut(l$x, val, labels=FALSE, include.lowest =TRUE, right=TRUE) 
    #xout <- sort(c(seq(val[1], val[length(val)], l = 150), val[-c(1, length(val))],val[-c(1, length(val))]*1.0001))
    #l <- suppressWarnings(approx(as.vector(m[, c(1, 3)]), as.vector(m[, c(2, 4)]), xout = xout))
    #val[length(val)]<-max(l$x) #aggiunto 11/09/17
    #id.group <- cut(l$x, val, FALSE, TRUE)
    yhat <- l$y
    xhat <- l$x
    m[, c(2, 4)] <- x$family$linkinv(m[, c(2, 4)])
    if (!add) {
      plot(as.vector(m[, c(1, 3)]), as.vector(m[, c(2,
                                                    4)]), type = "n", xlab = xlabs, ylab = ylabs,
           main = opz$main, sub = opz$sub, 
           cex.axis = opz$cex.axis,
           cex.lab = opz$cex.lab,
           xlim = opz$xlim,
           ylim = if(is.null(opz$ylim)) enl.range(fit, rangeCI, enlarge=dens.rug) else opz$ylim )
      if(dens.rug){
        density <- density(xvalues)
        # the height of the densityity curve
        max.density <- max(density$y)
        # Get the boundaries of the plot to
        # put the density polygon at the x-line
        plot_coordinates <- par("usr")
        # get the "length" and range of the y-axis
        y.scale <- plot_coordinates[4] - plot_coordinates[3]
        # transform the y-coordinates of the density
        # to the lower 10% of the plotting panel
        density$y <- (0.1 * y.scale / max.density) * density$y + plot_coordinates[3]
        ## plot the polygon
        polygon( density$x , density$y , border = FALSE , col = dens.col) 
        box()
      }
      
      if(rug) {
        #usare rug()?  
        segments(xvalues, rep(par()$usr[3],length(xvalues)), xvalues,
                 rep(par()$usr[3],length(xvalues))+ abs(diff(par()$usr[3:4]))/80)
      }
    }
    if (res) {
      if(hide.zeros) {
        fit <- fit[abs(xvalues)>1e-8]
        xvalues <- xvalues[abs(xvalues)>1e-8]
      }
      if(is.null(smoos)) { smoos <- if(length(xvalues)>10000) TRUE else FALSE }
      if(smoos){
        smoothScatter(xvalues, fit, add=TRUE, nrpoints = 0, colramp= colorRampPalette(c("white", res.col)))
      } else {
        points(xvalues, fit, cex = cexs, pch = pchs, col = res.col)
        
      }
    }
    if(conf.level>0){
      if(rev.sgn) vall<- -vall
      if(shade) {
        polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])),
                col = col.shade, border=NA) 
      } else {
        #browser()
        id.group1 <- cut(vall, val, labels=FALSE, include.lowest =TRUE, right=TRUE) #serve per gli IC..
        for (i in 1:max(id.group1)) matlines(vall[id.group1 == i], ciValues[id.group1 == i,-1], type="l", lty=2, col=cols[i])
        #matlines(vall, ciValues[,-1], type="l", lty=2, col=cols)
      }
    }
    
    yhat <- x$family$linkinv(yhat)
    if (length(cols) == 1) cols <- rep(cols, max(id.group))
    if (length(lwds) == 1) lwds <- rep(lwds, max(id.group))
    if (length(ltys) == 1) ltys <- rep(ltys, max(id.group))
    for (i in 1:max(id.group)) {
      lines(xhat[id.group == i], yhat[id.group == i], col = cols[i],
            lwd = lwds[i], lty = ltys[i])
      if(prev.trend) lines(xhat[xhat>est.psi[i]], x$family$linkinv((a[i]+b[i]*xhat)[xhat>est.psi[i]]), col=cols[i], lwd = lwds[i]*.65, lty = 2)
    }
    #-------------------------------------------------------------------------------
  } else { #se LM o "GLM con link=TRUE (ovvero linkinv=FALSE)"
    ##---> transf!!!
    y.val<- do.call(transf, list(y.val)) 
    y.val1<-do.call(transf, list(y.val1))
    r <- cbind(val, y.val)
    r1 <- cbind(val, y.val1)
    rr <- rbind(r, r1)
    fit <- c(y.val, y.val1)
    if (res) {
      ress <- if (inherits(x, what = "glm", which = FALSE))
        residuals(x, "working") #* sqrt(x$weights) mgcv::gam() usa " ..*sqrt(x$weights)/mean(sqrt(x$weights))"
      else resid(x)
      #if(!is.null(x$offset)) ress<- ress - x$offset
      #fit <- broken.line(x, term, gap = show.gap, link = link, interc = TRUE) + ress + const
      #fit <- predict.segmented(x, ifelse(rep(rev.sgn, length(xvalues)),-xvalues,xvalues), type=tipo) + ress + const
      fit <- fit0 + ress + const
    }
    if (!add)
      plot(rr, type = "n", xlab = xlabs, ylab = ylabs,
           main = opz$main, sub = opz$sub, 
           xlim = opz$xlim,
           cex.axis = opz$cex.axis,
           cex.lab = opz$cex.lab,
           #ylim = if(is.null(opz$ylim)) enl.range(fit, rangeCI, enlarge=dens.rug) else opz$ylim)
           ylim = if(is.null(opz$ylim)) enl.range(fit, rangeCI, do.call(transf, list(m[, c(2,4)])), enlarge=dens.rug) else opz$ylim)
    if(dens.rug){
      density <- density(xvalues)
      # the height of the densityity curve
      max.density <- max(density$y)
      # Get the boundaries of the plot to
      # put the density polygon at the x-line
      plot_coordinates <- par("usr")
      # get the "length" and range of the y-axis
      y.scale <- plot_coordinates[4] - plot_coordinates[3]
      # transform the y-coordinates of the density
      # to the lower 10% of the plotting panel
      density$y <- (0.1 * y.scale / max.density) * density$y + plot_coordinates[3]
      ## plot the polygon
      polygon(density$x , density$y , border = F , col = dens.col) 
      box()
    }
    if(rug) {segments(xvalues, rep(par()$usr[3],length(xvalues)), xvalues,
                      rep(par()$usr[3],length(xvalues))+ abs(diff(par()$usr[3:4]))/80)}
    if (res) {
      if(hide.zeros) {
        fit <- fit[abs(xvalues)>1e-8]
        xvalues <- xvalues[abs(xvalues)>1e-8]
      }
      if(is.null(smoos)) { smoos <- if(length(xvalues)>10000) TRUE else FALSE }
      if(smoos){
        smoothScatter(xvalues, fit, add=TRUE, nrpoints = 0, colramp= colorRampPalette(c("white", res.col)))
      } else {
        #browser()
        points(xvalues, fit, cex = cexs, pch = pchs, col = res.col)              
      }
    }
    if(rev.sgn) vall<- -vall
    if(conf.level>0) {
      if(shade) {
        polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])), col = col.shade, border=NA) 
      } else {
        #infittire vall, soprattutto in prossimita' dei psi?
        id.group1 <- cut(vall, val, labels=FALSE, include.lowest =TRUE, right=TRUE) #serve per gli IC..
        for (i in 1:max(id.group1)) matlines(vall[id.group1 == i], ciValues[id.group1 == i,-1], type="l", lty=2, col=cols[i])
        #VECCHIO: matlines(vall, ciValues[,-1], type="l", lty=2, col=cols)              
      }
    }
    #aggiunto 06/2019 perche' sotto disegnava linee (e non curve)
    #        segments(m[, 1], do.call(transf, list(m[, 2])), m[, 3], do.call(transf, list(m[, 4])), 
    #                col = cols, lwd = lwds, lty = ltys)
    #---
    # modificato 8/2/21.. adesso le linee si uniscono sempre.
    #.. con valori tipo 2010 (date), non si uniscono..
    #comunque vall ha piu' valori di xout, quindi e' sufficiente assegnare xout<-vall (01/10/2021)
    #xout <- sort(c(seq(val[1], val[length(val)], l = 50), val[-c(1, length(val))], 
    #              pmax(val[-c(1, length(val))]*1.0001, val[-c(1, length(val))]*.9999)))
    #if(rev.sgn) vall<- -vall
    xout <- vall
    l <- suppressWarnings(approx(as.vector(m[, c(1, 3)]), as.vector(m[, c(2, 4)]), xout = xout))
    val[length(val)]<- if(rev.sgn) min(l$x) else max(l$x) #aggiunto 11/09/17; messo il if .. else 9/3/21
    
    #id.group <- cut(l$x, val, labels=FALSE, include.lowest =TRUE, right=TRUE)
    id.group <- cut(vall, val, labels=FALSE, include.lowest =TRUE, right=TRUE) #e' come id.group1
    #---
    xhat <- l$x
    yhat <- l$y
    yhat <- do.call(transf, list(yhat)) #transf(yhat)
    if (length(cols) == 1) cols <- rep(cols, max(id.group))
    if (length(lwds) == 1) lwds <- rep(lwds, max(id.group))
    if (length(ltys) == 1) ltys <- rep(ltys, max(id.group))
    for (i in 1:max(id.group)) {
      lines(xhat[id.group == i], yhat[id.group == i], col = cols[i], lwd = lwds[i], lty = ltys[i])
      #if(conf.level>0 && !shade) matlines(vall[id.group1 == i], ciValues[id.group1 == i,-1], type="l", lty=2, col=cols[i])
      if(prev.trend) lines(xhat[xhat>est.psi[i]], (a[i]+b[i]*xhat)[xhat>est.psi[i]], col=cols[i], lwd = lwds[i]*.65, lty = 2)
    }
    #        if(prev.trend){
    #         for(i in 1:(length(est.psi)+1)) lines(xhat[xhat>est.psi[i]], a[i]+b[i]*xhat)[xhat>est.psi[i]], col=cols[i], lwd = lwds[i]*.7, lty = 2)
    #      }
  }
  invisible(NULL)
  }
}



