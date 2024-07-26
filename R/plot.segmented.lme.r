plot.segmented.lme<-function(x, level=1, id = NULL, res = TRUE, pop = FALSE,
         yscale = 1, xscale = 1, n.plot, pos.leg = "topright", vline = FALSE, lines = TRUE, 
         by=NULL, add=FALSE, conf.level=0, withI=TRUE, vcov.=NULL, shade=FALSE, drop.var=NULL, text.leg=NULL, 
         id.name=TRUE, ...) {
  # plotting fitted segmented relationships for multiple subjects
  # obj: a 'segmented.lme' object 
  # id: the subject id to be plotted
  # n.plot: a vector to be passed to par(mfrow= (should be coherent with length(id)). If missing, it is computed \epending on
  # length(id). 
  # If add=TRUE, it is expected that the plotting area is already split in different plots..  
  #         if prod(n.plot)<length(id), just the *first* `prod(n.plot)' units will be drawn. 
  # pop: if TRUE the population (fixed effects only) lines are added. In this case it could be useful to set the
  # same x-range (use xlim 
  # leg: where the id should be printed. Set to NULL for no id on the plot. yscale=1 => y range for all
  # subjects xscale=1 => x range for all subjects 
  # ...: arguments to pass to plotSegLme(), especially col.l, lwd.l, lty.l for the
  #     individual lines, and col.p, lty.p, lwd.p for the population
  #     lines (provided 'pop=T').  NB 'col' refers to points (provided that 'res=TRUE')
  #     col, col.l and col.p may be *vectors* (will be recycled)
  # lines: if FALSE, points (rather than lines) are plotted (useful if the segmented profile depends on
  #     additional covariates and cannot be displayed) x11() #quartz()?
  #vline:  if TRUE, a (dashed) vertical line is drawn to emphasize the individual breakpoint.
  #========================================================================
  return.psi <- function(obj, level) {
    # restituisce predizioni dei RE cahngepoints per diversi livelli
    # di nested. NON completa
    #--------------------
    fn.re <- function(obj) {
      # restituisce un array n x n.ranef x terms n ? il n. totale
      # delle misurazioni..  n.ranef ? il n. dei random effects
      # (tipicamente ? 1, >1 con nested..)  terms ? il n. dei
      # termini coinvolti nei random effects (ad es., intercept, x
      # ..)
      ro <- ranef(obj)
      n.levels <- ncol(obj$groups)  #n. dei livelli casuali (ad es., se nested..)
      if (n.levels <= 1) {
        ro <- list(ro)
        names(ro) <- names(obj$groups)
      }
      nomi.levels <- names(obj$groups)  #nomi degli effetti casuali names(ranef(obj))
      n.terms <- sapply(ro, ncol)
      nomiTermini <- unique(as.vector(unlist(sapply(ro, colnames))))
      tutti <- array(0, c(nrow(obj$groups), ncol(obj$groups), max(n.terms)),
                     dimnames = list(NULL, names(obj$groups), nomiTermini))
      for (nome in nomiTermini) {
        for (j in nomi.levels) {
          if (nome %in% names(ro[[j]])) {
            for (i in unique(obj$groups[, j])) tutti[obj$groups[,
                                                              j] == i, j, nome] <- ro[[j]][rownames(ro[[j]]) == i,
                                                                                             nome]
          }
        }
      }
      tutti
    }
    #---------------------
    if (missing(level)) level <- ncol(obj[[2]]$groups)
    psi <- rowSums(obj$misc$matrix.psi[, 1:(level + 1), drop = FALSE])
    if (level == 0) {
      names(psi) <- NULL
    } else {
      if (level < ncol(obj[[2]]$groups))
        names(psi) <- sapply(strsplit(names(psi), "/"), function(.x) .x[max(level,
                                                                            1)])
    }
    return(psi)
    # obj<-obj[[1]] RE<-fn.re(obj) if(! 'G0' %in% dimnames(RE)[[3]])
    # stop('no 'G0' term!') n.ranef<-dim(RE)[2] psi.i<-vector('list',
    # ncol(obj$groups)) for(j in 1:n.ranef) { valori.sum<-
    # rowSums(RE[, 1:j,'G0',drop=FALSE]) valori<- tapply(valori.sum,
    # obj$groups[,j], function(x)x[1])
    # names(valori)<-names(tapply(obj$groups[,j], obj$groups[,j],
    # function(x)x[1])) psi.i[[j]]<-valori } names(psi.i) <-
    # colnames(obj$groups) attr(psi.i, 'grpNames') <-
    # attr(ranef(obj), 'grpNames') return(psi.i) kappa0i <- kappa0+ki
    # etai<-kappa0i if(id.z.psi) { kappa.old<-kappa #length=1
    # kappa<-fixed.effects(obj)[nomiG] #esclude G0..
    # etai<-etai+drop(Z.psi%*%kappa) } psi.ex<-if(psi.link=='logit')
    # inv.logit(etai,min.Z,max.Z) else etai
  }
  #========================================================================
  plotSegLme <- function(obj, id = stop("'id' should be provided"), add = FALSE, res = TRUE,
                         pop = FALSE, yscale = -1, xscale = -1, text.leg = paste("id =", id), pos.leg = NULL, vline = FALSE,
                         xLab, yLab, level, lines = TRUE, opzione = 1, ...) { #line.col=1, res.col=grey(.7),
    # Simply plots (or adds) the observed data and the segmented fitted lines for subject
    # 'id'
    #---
    # obj: an object of class 'segmented.lme' id: the subject 'id' add: if FALSE, a new
    # plot is produced with observations and fitted lines superimposed.  res: if TRUE the
    # observations (partial residuals) are added; otherwise only the fitted lines pop: if
    # TRUE the population-level estimate of the segmented relationship is added.. yscale
    # if <0, the y-scale refers to the values of 'id' only; otherwise the overall range
    # relevant to *all* subjects (useful for comparisons) main: the plot title. It can be
    # '' leg: if !NULL it can be one of 'top', 'topright',... and the id subject is put
    # on the plot. vline: if TRUE lines: if FALSE, points (rather than lines) are plotted
    # (useful if the segmented profile depends on additional covariates and cannot be
    # displayed) 
    #...: argomenti da passare al plot, compresi 'col.l' e 'lwd.l' 'lty.l'
    # che servono per le segmented lines individuali e col.p, lty.p, lwd.p che servono
    # per le linee di pop (se pop=TRUE) Problema: se ci sono nested re i levels dell
    # innermost factor vengono modificati 'lev1/lev2' (che rappresentano i nomi dei
    # coef() o ranef) dove lev1 e' il livello dell'altro fattore. Quindi o in 'id' si
    # specifica il livello costruito 'lev1/lev2' oppure i nomi dei coef() si devono
    # modificare estraendo solo il valore dopo il '\'
    rnfGrp <- obj$lme.fit.noG$groups
    if (missing(level)) level <- ncol(rnfGrp)
    # con nested funziona ma i valori dei fitted sono 'strani..', cio? dovrebbero essere
    # comunque una relazione segmented S? ? giusto solo i psi sono diversi...
    if (level != ncol(rnfGrp))
      stop("Currently only the innermost level (i.e., ", ncol(rnfGrp), ") is allowed")
    if (level > ncol(rnfGrp) || level < 0)
      stop(" 'level' should be an integer in [0, ", ncol(rnfGrp), "]")
    # rownames(coef(obj[[2]])) sono the innermost levels
    nomi <- levels(rnfGrp[, max(1, level)])  #level 0 non funziona
    if (!(id %in% nomi))
      stop("the specified 'id' is not consistent with 'level'  (should be ", nomi[1], ",",
           nomi[2], ",..,", nomi[length(nomi)], ")")
    # psi<- obj$psi.i[[paste(id)]] in realt? si possono ottenere i psi a diversi livelli
    # di nested
    psi <- if (ncol(obj$misc$matrix.psi) <= 1)
      obj$misc$matrix.psi[, 1] else rowSums(obj$misc$matrix.psi[, 1:(level + 1)])
    # se ci sono nested, i nomi di psi saranno sempre con '/' che non e' coerente se
    # level=1
    if (level < ncol(rnfGrp)) {
      if (ncol(rnfGrp) >= 2) {
        names(psi) <- sapply(strsplit(names(psi), "/"), function(.x) .x[level])
      }
    }
    psi <- psi[[paste(id)]]
    nameID <- names(rnfGrp)[level]
    y <- resid(obj$lme.fit.noG) + fitted(obj$lme.fit.noG)  #
    if (level < ncol(rnfGrp))
      names(y) <- obj$lme.fit.noG$groups[, max(level, 0)]
    range.ok <- range(y)
    XY <- cbind(obj$Z, y)
    x <- XY[names(y) == id, 1]
    y <- XY[names(y) == id, 2]
    if (yscale < 0) range.ok <- range(y)
    
    #browser()
    
    range.ok[1] <- if (sign(range.ok[1]) > 0)
      range.ok[1] * 0.98 else range.ok[1] * 1.03
    range.ok[2] <- if (sign(range.ok[2]) > 0)
      range.ok[2] * 1.03 else range.ok[2] * 0.98
    
    # x<-obj$Z
    rangeX.ok <- range(obj$Z)
    # x<-x[names(obj$Z)==id]
    if (xscale < 0) rangeX.ok <- range(x)
    
    opz <- list(...)
    
    #browser()
    
    
    opz$x <- x
    opz$y <- y
    
    #browser()
    
    if (is.null(opz$col)) opz$col <- grey(.7) 
    if (is.null(opz$cex)) opz$cex <- 1.5
    if (is.null(opz$pch)) opz$pch <- 19
    if (is.null(opz$ylim)) opz$ylim <- range.ok
    if (is.null(opz$xlim)) opz$xlim <- rangeX.ok
    if (!res) opz$type <- "n"
    
    
    if (is.null(opz$ylab)) opz$ylab<- all.vars(formula(obj[[2]]))[1]
    if (is.null(opz$xlab)) opz$xlab<- obj$namesGZ$nameZ
    
    
    
    #if (missing(yLab)) yLab <- "response"
    #if (missing(xLab)) xLab <- obj$namesGZ$nameZ
    #opz$ylab <- yLab
    #opz$xlab <- xLab
    
    
    # set col.p, lwd.p, lty.p for *population* lines (provided pop=TRUE)
    if (!is.null(opz$p.col)) {
      p.col <- opz$p.col
      opz$p.col <- NULL
    } else {
      p.col <- 1
    }
    
    if (!is.null(opz$p.lty)) {
      lty.p <- opz$p.lty
      opz$p.lty <- NULL
    } else {
      p.lty <- 2
    }
    
    if (!is.null(opz$p.lwd)) {
      p.lwd <- opz$p.lwd
      opz$p.lwd <- NULL
    } else {
      p.lwd <- 1.5
    }
    
    # set col.l, lwd.l, lty.l for *individual* lines (or points if 'lines=FALSE')
    if (!is.null(opz$l.pch)) {
      l.pch <- opz$l.pch
      opz$l.pch <- NULL
    } else {
      l.pch <- 3
    }
    
    if (!is.null(opz$l.lty)) {
      l.lty <- opz$l.lty
      opz$l.lty <- NULL
    } else {
      l.lty <- 1
    }
    
    if (!is.null(opz$l.col)) {
      l.col <- opz$l.col
      opz$l.col <- NULL
    } else {
      l.col <- 1
    }
    if (!is.null(opz$l.lwd)) {
      l.lwd <- opz$l.lwd
      opz$l.lwd <- NULL
    } else {
      l.lwd <- 2
    }
    
    if (!is.null(opz$t.col)) {
      t.col <- opz$t.col
      opz$t.col <- NULL
    } else {
      t.col <- 1
    }
    
    
    
    
    if (!add) do.call(plot, opz)
    if (!is.null(pos.leg)){ 
      legg <- if(id.name) paste(nameID, id, sep = " = ") else id
      legend(pos.leg, legend = legg, bty = "n", text.col=t.col)
    }
    
    
    ff <- fitted(obj$lme.fit.noG, level = level)
    mu <- ff[names(ff) == id]  #? fitted.segmented.lme(fit,1)
    
    m <- cbind(x, mu)
    m <- m[order(m[, 1]), ]
    
    
    #browser()
    
    
    if (!lines) {
      points(m, col = l.col, pch = l.pch, lwd = l.lwd)  #do.call(points, opz)
    } else {
      if (opzione == 1) {
        mL <- m[m[, 1] <= psi, , drop = FALSE]
        if (nrow(mL) > 1) {
          # if(length(unique(round(diff(mL[,1])/diff(mL[,2]), 6)))>1)
          # stop('Apparently non-segmented profile for unit #', id, '.. try
          # 'lines=FALSE'\n')
          fL <- splinefun(mL[, 1], mL[, 2])
          f.psi <- fL(psi)
        } else {
          mR <- m[m[, 1] >= psi, , drop = FALSE]
          # if(length(unique(round(diff(mR[,1])/diff(mR[,2]), 6)))>1)
          # stop('Apparently non-segmented profile for unit #', id,'.. try
          # 'lines=FALSE'\n')
          fR <- splinefun(mR[, 1], mR[, 2])
          f.psi <- fR(psi)
        }
        # lines(c( mL[,1], psi, mR[,1]), c( mL[,2], f.psi, mR[,2]), col=col.l,
        # lwd=lwd.l, lty=lty.l)
        
        lines(c(m[1, 1], psi, m[nrow(m), 1]), c(m[1, 2], f.psi, m[nrow(m), 2]), col = l.col,
              lwd = l.lwd, lty = l.lty)
      } else {
        nomiCoef.ok <- intersect(c("(Intercept)", obj$namesGZ$nameZ, "U"), colnames(coef(obj$lme.fit.noG)))
        nomiCoef.ok <- c(nomiCoef.ok, obj$namesGZ$nomiUx)
        coef.ok <- coef(obj$lme.fit.noG, level = level)[id, nomiCoef.ok]
        psi.ok <- return.psi(obj, level = level)[id]  #dipende da eventuali termini fissi in z.psi, e dal livello di nesting.
        # adesso devi copstruire la matrice del disegno. 
        
        
        #########################################################
        ### GUARDARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #########################################################
        
        #browser()
        
        xvar<-seq(min(x), max(x), l=100)
        
        coef.ok<- as.numeric(coef.ok)
        
        if(!is.null(obj$namesGZ$nomiUx)) {
          coef.ok[3] <- coef.ok[3]+ coef.ok["bUx"]* obj$lme.fit.noG$data[,obj$namesGZ$nomiUx]
        }
        
        mu.ok <- cbind(1, xvar, pmax(xvar-psi.ok,0))%*%coef.ok
        lines(xvar, mu.ok, col = l.col, lwd = l.lwd, lty = l.lty)
        
        
        # come fare se ci sono variabili U.x? coef.ok include il coef ma c'? bisogno
        # del valore corrispondente all'unita' id..
        # obj$lme.fit.noG$data[,obj$namesGZ$nomiUx] e poi selezionare per l'unita'
        # 'id' E se ci sono interazioni con intercetta e left slope???
      }
    }
    if(vline) segments(psi, par()$usr[3], psi, f.psi, lty = 3, col = l.col)
    
    #browser()
    if(attr(obj$psi.i,"is.break")[paste(id)]){
      points(psi, par()$usr[3] * 1, pch = "X", col = l.col, cex = 1.2)
      points(psi, f.psi, pch = "x", col = l.col, cex = 1, lwd=1.5)
    }
    #browser()
    
    # codici vecchi..  #left side mL<-m[m[,1]<=psi, ,drop=FALSE] fL<-splinefun(mL[,1],
    # mL[,2]) new.xL<- c(min(mL[,1]), psi) #right side mR<-m[m[,1]>=psi, ,drop=FALSE]
    # fR<-splinefun(mR[,1], mR[,2]) new.xR<- c(psi, max(mR[,1])) lines(new.xL,
    # fL(new.xL), col=1, lwd=2) lines(new.xR, fR(new.xR), col=1, lwd=2) if(vline)
    # segments(psi, par()$usr[3], psi, fR(psi), lty=3, col=1)
    if (pop) {
      #browser()
      # mu<-fitted(obj[[2]])[names(fitted(obj[[2]]))==id] #e'
      # fitted.segmented.lme(fit,1)
      # mu<-fitted(obj[[2]],0)[names(fitted(obj[[2]],0))==id] #e'
      # fitted.segmented.lme(fit,0)
      
      # mu<-fitted(obj, level=0)[names(fitted(obj, level=0))==id] #funziona solo se i
      # dati sono ordinati (osservazioni dello stesso individuo vicine..). Per cui il
      # 14/7 messo la seguente: browser() mu<-fitted(obj[[2]],
      # level=0)[names(fitted(obj[[2]], level=ncol(obj[[2]]$groups)))==id]
      mu <- fitted(obj, level = 0)
      mu <- mu[names(mu) == id]
      # mu<-mu[names(fitted(obj, level=ncol(obj[[2]]$groups)))==id]
      psi <- obj$fixed.psi[[paste(id)]]
      m <- cbind(x, mu)
      m <- m[order(m[, 1]), ]
      
      # mL<-m[m[,1]<=psi, ,drop=FALSE] if(nrow(mL)>1){ fL<-splinefun(mL[,1], mL[,2])
      # f.psi<-fL(psi) } else { mR<-m[m[,1]>=psi, ,drop=FALSE] fR<-splinefun(mR[,1],
      # mR[,2]) f.psi<-fR(psi) } lines(c( m[1,1], psi, m[nrow(m),1]), c( m[1,2], f.psi,
      # m[nrow(m),2]), col=col.l, lwd=lwd.l) left side
      m1 <- m[m[, 1] <= psi, , drop = FALSE]
      # right side
      m2 <- m[m[, 1] >= psi, , drop = FALSE]
      if (nrow(m1) > 0) {
        f1 <- splinefun(m1[, 1], m1[, 2])
        estremo <- if (nrow(m2) > 0)
          psi else min(psi, max(m1[, 1]))
        new.x1 <- c(min(m1[, 1]), estremo)
      }
      
      if (nrow(m2) > 0) {
        f2 <- splinefun(m2[, 1], m2[, 2])
        # new.x1<- seq(psi, max(m1[,1]), l=200)
        estremo <- if (nrow(m1) > 0)
          psi else max(psi, min(m2[, 1]))
        new.x2 <- c(estremo, max(m2[, 1]))
      }
      if (nrow(m1) > 0) {
        if (nrow(m1) > 1)
          lines(new.x1, f1(new.x1), col = p.col, lwd = p.lwd, lty = p.lty) else lines(new.x1, c(f1(new.x1)[1], f2(new.x2)[1]), col = p.col, lwd = p.lwd, lty = p.lty)
      }
      if (nrow(m2) > 0) {
        if (nrow(m2) > 1)
          lines(new.x2, f2(new.x2), col = p.col, lwd = p.lwd, lty = p.lty) else lines(new.x2, c(f1(new.x1)[2], f2(new.x2)[2]), col = p.col, lwd = p.lwd, lty = p.lty)
      }
      points(psi, par()$usr[3] * 1.015, pch = 4, col = p.col)
      # segments(psi, par()$usr[3], psi, f1(psi), lty=3, col=1)
    }
  }
  #========================================================================
  plotmarg<-function(obj, by=NULL, add=FALSE, conf.level=0, pos.leg=NULL, withI=TRUE, vcov.=NULL, 
                     shade=FALSE, drop.var=NULL, text.leg, ... ){
    #=========>> da provare con by con piu' termini e se leftSlope=0 
    #obj: the segmented.lme object
    #by: a named list indicating covariate names and corresponding values affecting the fitted segmented relationship.
    #    Example: by=list(group="2",z2=.2) 
    #conf.level: for pointwise CI..
    #pos.leg if different from NULL, a legend is added on the plot at "pos.leg"
    #  (which should be one a single keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center")
    #withI: if TRUE, the fitted lines are plotted with intercept (if included in the model)
    #vcov.: the fixed effect cov matrix. If NULL is computed by vcov.segmented.lme
    #shade: logical (ignored if conf.level=0)
    #drop.var: possible coefficient names to be removed before computing the segmented relationship (E.g. the group-specific intercept..)
    #...: further arguments (col, lty, lwd, xlim, ylim, xlab,..) to be passed to matplot()/matlines(). Can be vectors (e.g. lty=c(2,1,2) and lwd=c(1,2,1)
    #text.leg: if specified, it is legend to be added, provided pos.leg has been specified. 
    
    #obj[[1]] -> obj$lme.fit
    #obj[[2]] -> obj$lme.fit.noG
    #------------------
    if((conf.level<0 || conf.level>=1)) stop(" 'conf.level' is meaningless'")
    if(conf.level>0){
      V<- if(is.null(vcov.)) vcov.segmented.lme(obj) else vcov. #object$lme.fit$varFix 
    } #s
    #  if(conf.level!=0 && !is.null(by)) stop("Not yet implemented..")
    if(!is.null(by) && !is.list(by)) stop("if provided, 'by' should be a (named) list of scalars")
    opz<-list(...)
    if("pos.legend" %in% names(opz)) warning("'pos.legend' ignored.. Do you mean 'pos.leg'?", call.=FALSE)
    if(!is.null(pos.leg)) pos.leg<- match.arg(pos.leg, c("bottomright", "bottom", "bottomleft", 
                                                         "left", "topleft", "top", "topright", "right", "center"))
    Z<-obj$Z 
    nomeZ<-obj$namesGZ$nameZ
    beta.noG<- fixef(obj$lme.fit.noG) 
    beta.all<-fixef(obj$lme.fit)
    beta.G<-beta.all[setdiff(names(fixef(obj$lme.fit)), names(beta.noG))]
    nomiCoef<-names(beta.noG)
    
    if(!is.null(by)) {
      a<-by
      #isZero<-sapply(a, function(x) x==0)
      if(!all(sapply(a, length)==1)) stop("vectors in 'by' are not allowed")
      nomiOK<-const<-idList<-vector("list", length(a))
      values<-vector(,length(a))
      
      for(i in 1:length(a)) {
        nomiOK[i]<-nomeOK <- if(is.character(a[[i]])) paste(names(a[i]),a[[i]], sep="") else names(a[i])
        
        #replace 0
        if(a[[i]]==0) a[[i]]<- 1e-16
        
        #per l'intercetta
        bInterc<- beta.noG[nomeOK]
        if(!is.character(a[[i]])) bInterc<-bInterc*a[[i]] 
        
        #per la left slope
        bLeftSlope<-c(beta.noG[paste(obj$namesGZ$nameZ,":",nomeOK, sep="")],
                      beta.noG[paste(nomeOK, ":", obj$namesGZ$nameZ, sep="")])
        bLeftSlope<-bLeftSlope[!is.na(bLeftSlope)]
        if(!is.character(a[[i]])) bLeftSlope<-bLeftSlope*a[[i]]
        if(length(bLeftSlope)<=0) bLeftSlope<-NA 
        
        #per la slope-diff
        bU<-beta.noG[paste("U", nomeOK, sep=".")]
        if(!is.character(a[[i]])) bU <-bU*a[[i]]
        
        #per il changepoint
        bG<-beta.G[paste("G", nomeOK, sep=".")]
        if(!is.character(a[[i]])) bG <-bG*a[[i]]
        
        const[[i]]<-c(bInterc, bLeftSlope, bU, bG)
        const[[i]]<- ifelse(is.na(const[[i]]),0,const[[i]])
        
        idList[[i]]<-names(c(bInterc, bLeftSlope, bU, bG))
        values[i]<-ifelse(is.character(a[[i]]),1,a[[i]])
      }
      #browser()          
      const<-matrix(unlist(const),4, byrow=FALSE)
      colnames(const)<-names(by)
      nomiNOdiff <- names(which(colSums(const)==0))
      if(length(nomiNOdiff)>0) warning("The", paste(" '", nomiNOdiff,"' ",sep=""), "value supplied in 'by' does not modify the baseline line", call. = FALSE)
      
      nomiCoef<- c("(Intercept)", nomeZ, "U", "G0", unlist(idList))
      ##########################################    
    } else { #se 'by' e' NULL
      const<-matrix(0,4,1) 
      nomiCoef<- c("(Intercept)", nomeZ, "U", "G0")
      values<-rep(1,4)
    }  
    ##########################################    
    #browser()  
    #prepara la matrice del disegno..
    est.psi.fixed <-fixef(obj$lme.fit)["G0"]+ sum(const[4,])
    if(obj$call$psi.link=="logit") est.psi.fixed <- plogis(est.psi.fixed)
    Z.new<-as.numeric(sort(c(seq(min(Z),max(Z),l=100), est.psi.fixed)))
    U<-pmax(Z.new-est.psi.fixed,0)
    X<-cbind(1,Z.new,U) #colnames(X)<-c("(Intercept)",nomeZ,"U")
    
    Ident<-diag(ncol(X))
    M<-vector("list", length=ncol(const))
    for(j in 1:ncol(const)) M[[j]]<-values[j]*Ident[, which(const[-4,j]!=0), drop=FALSE]
    
    M<-cbind(Ident, do.call("cbind", M))
    
    if(!withI) {
      M<-M[,-1] 
      nomiCoef <-setdiff(nomiCoef, "(Intercept)")
    }
    
    final.names<-setdiff(nomiCoef, c("G0",obj$namesGZ$nomiG,""))
    final.names<-final.names[!is.na(final.names)]
    
      #browser()
    
    if(!is.null(drop.var)){
      colnames(M)<-final.names
      final.names <-setdiff(final.names, drop.var)
      M<-M[, final.names]
    }
    
    XX<- X%*%M
    r<-fit<- XX %*% beta.noG[final.names]
    
    if (conf.level > 0) {
      zalpha<- -qnorm((1-conf.level)/2)
      V<-V[final.names,final.names]
      SE.fit<-sqrt(rowSums((XX %*% V) * XX)) #sqrt(diag(X%*%Var%*%t(X)))
      r<-cbind(fit-zalpha*SE.fit, fit, fit+zalpha*SE.fit)
      opz$lty<-c(2,1,2)
    }
    
    #  b<-beta.noG[c("(Intercept)",nomeZ,"U")]
    #  b<-b+ rowSums(const[1:3,,drop=FALSE])
    #  b<-b[!is.na(b)]
    #  r<-fit<- drop(X  %*% b)
    if(is.null(opz$col))  opz$col<-"1"
    if(is.null(opz$type)) opz$type<-"l"
    if(is.null(opz$lty))  opz$lty<-1
    if(is.null(opz$lwd))  opz$lwd<-1.8
    if(is.null(opz$xlab)) opz$xlab<-obj$namesGZ$nameZ
    if(is.null(opz$ylab)) opz$ylab<-all.vars(formula(obj[[1]]))[1]
    if(!is.null(opz$alpha.f)) {alpha.f<-opz$alpha.f;opz$alpha.f<-NULL} else {alpha.f<-.15}  
    opz$x<-Z.new
    opz$y<-r
    if(add) do.call(matlines, opz)  else do.call(matplot, opz)
    
    if (shade && conf.level>0) polygon(c(Z.new, rev(Z.new)), c(r[, 1], rev(r[, 3])), 
                                       col = adjustcolor(opz$col, alpha.f), border = NA)
    
    #browser()  
    if(!is.null(pos.leg)) {
      if(!is.null(by)){
        id<-apply(const,2,function(.x) any(.x!=0)) #funziona sempre?
        #se sono tutti 0, significa che stai disegnando la segmented baseline.. 
        #Ma perche' questo dovrebbe influenzare la scelta di mettere la legenda? proviamo a metterlo sempre TRUE
        #id<-TRUE
        leg<-paste(names(by)[id],unlist(by)[id],sep="=",collapse="  ")
      } else {
        leg<-NULL
      }
      if(!is.null(text.leg)) leg <- text.leg
      try(legend(pos.leg, leg, bty="n", col=min(opz$col), lty=min(opz$lty), lwd=max(opz$lwd)), silent=TRUE)
    }
  }
  
  #========================================================================
  #inizio funzione
  if(level==0){
    plotmarg(x, by=by, add=add, conf.level=conf.level, pos.leg = pos.leg, withI=withI, vcov.=vcov., 
             shade=shade, drop.var=drop.var, text.leg=text.leg, ...)
  } else {
    obj <- x
    rnfGrp <- obj$lme.fit.noG$groups
    if (missing(level))
      level <- ncol(rnfGrp)
    if (level != ncol(rnfGrp))
      stop("Currently only the innermost level (i.e., ", ncol(rnfGrp), ") is allowed")
    opz <- list(...)
    Ylab <- if (is.null(opz$ylab))
      all.vars(formula(obj[[1]]))[1] else opz$ylab
    Xlab <- if (is.null(opz$xlab))
      obj$namesGZ$nameZ else opz$xlab
    opz$xlab <- opz$ylab <- NULL
    opz$pop <- pop
    opz$res <- res
    opz$xLab <- ""
    opz$yLab <- ""
    opz$main <- ""
    opz$xaxt <- "n"
    opz$yaxt <- "n"
    opz$pos.leg <- pos.leg
    opz$yscale <- yscale
    opz$xscale <- xscale
    opz$vline <- vline
    opz$level <- level
    opz$lines <- lines
    opz$obj <- quote(obj)
    if (is.null(id))
      id <- names(obj$psi.i)  #levels(rnfGrp[,ncol(rnfGrp)])  
    if (missing(n.plot))
      n.plot <- if (length(id) <= 1)
        c(1, 1) else c(3, ceiling(length(id)/3))
    
    if(prod(n.plot)!=1) id <- id[1:min(prod(n.plot), length(id))]
    
    # color of individual lines
    l.col <- if (!is.null(opz$l.col)) opz$l.col else 1
    l.col <- rep(l.col, length(id))
    
    l.lwd <- if (!is.null(opz$l.lwd)) opz$l.lwd else 1
    l.lwd <- rep(l.lwd, length(id))

    l.lty <- if (!is.null(opz$l.lty)) opz$l.lty else 1
    l.lty <- rep(l.lty, length(id))
    
    col <- if (!is.null(opz$col)) opz$col else grey(.7) #for residuals..
    col <- rep(col, length(id))
    
    p.col <- if (!is.null(opz$p.col)) opz$p.col else 1
    p.lty <- if (!is.null(opz$p.lty)) opz$p.lty else 3
    p.lwd <- if (!is.null(opz$p.lwd)) opz$p.lwd else 1
    #p.col <- rep(p.col, length(id))
    
    t.col <- if (!is.null(opz$t.col)) opz$t.col else 1 #for legend text..
    t.col <- rep(t.col, length(id))
    
     #browser()   
    
    # if(dev.cur()==1) { #se non e' aperto alcun device..
    #unico grafico con tutti i profili individuali.. 
    #======================================================
    xlim<- if(!is.null(opz$xlim)) opz$xlim else range(x$Z)
    ylim<- if(!is.null(opz$ylim)) opz$ylim else NULL
    k<-1
    if(prod(n.plot)==1 && length(id)>1){
      plotSegLme(obj=x, id=id[k], add=FALSE, pop=FALSE, res=FALSE, xLab='', yLab='', 
                 l.col=l.col[k], l.lwd=l.lwd[k], l.lty=l.lty[k], xlim=xlim, ylim=ylim)
      for (i in id[-1]) {
        k<-k+1
        plotSegLme(obj=x, id=i, add=TRUE, pop=FALSE, res=FALSE, xLab='', yLab='', 
                   l.col=l.col[k], l.lwd=l.lwd[k], l.lty=l.lty[k])
        # main='', xaxt='n', yaxt='n', leg=leg, yscale=yscale,
        # vline=vline, xscale=xscale, level=level, ...)
        #opz$id <- i
        #if (col.l.id) opz$col.l <- col.l[k]
        #if (col.p.id) opz$col.p <- col.p[k]
        #if (col.id) opz$col <- col[k]
        
        #guarda bene il discorso della stima sella relazione only-fixed-effects..
        #opz$pop<-FALSE
        #do.call(plotSegLme, opz)
      }
      box()
      if(pop) plotmarg(obj, add=TRUE, col=p.col, lwd=p.lwd, lty=p.lty)
      return(invisible(NULL))
    }
    
    old.mar<- par()$mar
    old.oma<- par()$oma
    old.mfrow<- par()$mfrow
    
    #======================================================
    if(length(id)>1){
        par(mfrow = n.plot)
        id.sx <- 1 + n.plot[2] * (0:(n.plot[1] - 1))  #i grafici di sx
        id.bot <- (prod(n.plot):1)[1:n.plot[2]]  #i grafici di sotto
        if(yscale>=0) {
          if(xscale>=0) {
            par(mar = rep(0, 4)) 
            } else {
              opz$xaxt<-"n"
              par(mar = c(0,0,2.2,0))
              }
        } else {
          #opz$yaxt<-NULL
          if(xscale>=0) {
            par(mar = c(0,0,0,2.5)) 
            } else {
              opz$xaxt<-"n"
              par(mar = c(0,0,2,2.5))
            }
        }
        par(oma = c(5, 5, 1, 1))
        out <- TRUE
      } else {
        id.sx <- 1:length(id)
        id.bot <- 1:length(id)
        out <- FALSE
      }
    k <- 0
    opz$add<- add
    #col.l lo prende anche sui punti????
    for (i in id) {
      #browser()
      k <- k + 1
      # plotSegLme(obj, id=i, pop=pop, res=res, xLab='', yLab='',
      # main='', xaxt='n', yaxt='n', leg=leg, yscale=yscale,
      # vline=vline, xscale=xscale, level=level, ...)
      opz$id <- i
      opz$l.col <- l.col[k]
      opz$p.col <- p.col[k]
      opz$col <- col[k]
      
      #guarda bene il discorso della stima sella relazione only-fixed-effects..
      opz$pop<-FALSE
      do.call(plotSegLme, opz)
      if(pop) plotmarg(obj, add=TRUE, lty=2)
      #browser()
      
      # tt<-axTicks(1) las=2
      #if((xscale>=0)&&(k %in% id.bot)) axis(1, cex.axis = 0.7, at = NULL) else axis(1, labels = FALSE)
      #if((yscale>=0)&&(k%in%id.sx)) axis(2, labels = TRUE, cex.axis = 0.7) else axis(2, labels = FALSE)
      #if((k %in% id.bot)) axis(1, cex.axis = 0.7, at = NULL) else axis(1, labels = FALSE)
      #
      #if(xscale>=0) {
       # if((k %in% id.bot)) axis(1, cex.axis = 0.7, at = NULL)
      #  axis(1, labels = FALSE)
      #}
      #browser()
      if((xscale<0)|| ((xscale>=0)&&(k%in%id.bot))) axis(1, labels = TRUE, cex.axis = 0.7) else axis(1, labels = FALSE)
      #if(((xscale>=0)&&(k%in%id.bot))) axis(1, labels = TRUE, cex.axis = 0.7) else axis(1, labels = FALSE)
      
      if((yscale<0)|| ((yscale>=0)&&(k%in%id.sx))) axis(2, labels = TRUE, cex.axis = 0.7) else axis(2, labels = FALSE)
      
    }
    if(length(id)>1) {
      mtext(Xlab, 1, line = 3, outer = out)
      mtext(Ylab, 2, line = 3, outer = out)
      par(mar=old.mar, oma=old.oma, mfrow=old.mfrow)
    }
  }
}
