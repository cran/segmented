`pscore.test`<- function(obj, seg.Z, k = 10, alternative = c("two.sided", "less", "greater"),
                         values=NULL, dispersion=NULL, df.t=NULL, more.break=FALSE, n.break=1, only.term=FALSE) { 
  #-------------------------------------------------------------------------------
  test.Sc2<-function(y, z, xreg, sigma=NULL, values=NULL, fn="pmax(x-p,0)", df.t="Inf", alternative, w=NULL, offs=NULL, 
                     nbreaks=1, ties.ok=FALSE, only.term=FALSE){
    #xreg: la matrice del disegno del modello nullo. Se mancante viene assunta solo l'intercetta.
    #Attenzione che se invXtX e xx vengono entrambe fornite, non viene fatto alcun controllo
    #invXtX: {X'X}^{-1}. if missing it is computed from xreg
    #sigma: the sd. If missing it is computed from data (under the *null* model)
    #values: the values with respect to ones to compute the average term. If NULL 10 values from min(z) to max(z) are taken.
    if(!is.null(offs)) y<-y-offs
    n<-length(y)
    if(missing(xreg)) xreg<-cbind(rep(1,n))
    id.ok<-complete.cases(cbind(y,z,xreg))
    y<-y[id.ok]
    z<-z[id.ok]
    xreg<-xreg[id.ok,,drop=FALSE]
    n<-length(y)
    k=ncol(xreg) #per un modello ~1+x
    if(is.null(values)) values<-seq(min(z), max(z), length=10)
    n1<-length(values)
    PSI<-matrix(values, nrow=n, ncol=n1, byrow=TRUE) #(era X2) matrice di valori di psi
    if(is.matrix(z)) {
      X1<-matrix(z[,1], nrow=n, ncol=n1, byrow=FALSE)
      X2<-matrix(z[,2], nrow=n, ncol=n1, byrow=FALSE)
      X<-eval(parse(text=fn), list(x=X1, y=X2, p=PSI)) #X<-pmax(X1-X2,0)
      pmaxMedio<-rowMeans(X)
    } else {
      X1<-matrix(z, nrow=n, ncol=n1, byrow=FALSE) #matrice della variabile Z
      if(length(fn)<=1){
        X<-eval(parse(text=fn), list(x=X1, p=PSI)) #X<-pmax(X1-PSI,0)
        pmaxMedio <- rowMeans(X)
        if(nbreaks>1){
          XX<-sapply(1:length(values), function(.x) X[,-(1:.x), drop=FALSE])
          XX<-do.call("cbind", XX)
          if(ties.ok) XX<-cbind(X, XX)
          pmaxMedio2 <- rowMeans(XX)
          pmaxMedio <- cbind(pmaxMedio, pmaxMedio2)
        }
      } else {
        pmaxMedio<-matrix(NA,n,length(fn))
        #list.X<-vector("list", length=length(fn))
        for(j in 1:length(fn)){
          #list.X[[j]]<-eval(parse(text=fn[j]), list(x=X1, p=PSI))
          X<-eval(parse(text=fn[[j]]), list(x=X1, p=PSI))
          pmaxMedio[,j]<-rowMeans(X)
        }
      }
    }
    if(only.term) return(pmaxMedio)
    if(is.null(w)) w<-1
    invXtX<-solve(crossprod(sqrt(w)*xreg))
    IA<-diag(n) - xreg%*%tcrossprod(invXtX, xreg*w) #I-hat matrix
    sc<-t(pmaxMedio*w) %*% IA  %*% y
    v.s<- t(pmaxMedio*w) %*% crossprod(t(IA)/sqrt(w))%*%(w*pmaxMedio)
    ris<-if(nbreaks==1) drop(sc/(sigma*sqrt(v.s))) else drop(crossprod(sc,solve(v.s,sc)))/(sigma^2)
    #if(length(fn)<=1 && cadj) ris<- sign(ris)*sqrt((ris^2)*(1-(3-(ris^2))/(2*n)))
    #passa alla F..
    df.t<-eval(parse(text=df.t))
    p2<- if(nbreaks==1) 2*pt(abs(ris), df=df.t, lower.tail=FALSE) else pchisq(ris, df=nbreaks, lower.tail=FALSE)#pf((ris/nbreaks)/(sigma^2), df1=nbreaks, df2=df.t, lower.tail =FALSE)#
    pvalue<-switch(alternative,
                less = pt(ris, df=df.t, lower.tail =TRUE) ,
                greater = pt(ris, df=df.t, lower.tail =FALSE) ,
                two.sided = p2)
    #pvalue<- 2*pt(abs(ris), df=df.t, lower.tail =FALSE) 
    r<-c(ris, pvalue)#, pmaxMedio)
    r
  #return(pmaxMedio)
    }
  
  #-------------------------------------------------------------------------------
  scGLM<-function(y, z, xreg, family, values = NULL, size=1, weights.var,
                  fn="pmax(x-p,0)", alternative=alternative){
    #score test for GLM
    #size (only if family=binomial())
    #weights.var: weights to be used for variance computations. If missing the weights come from the null fit
    output<-match.arg(output)
    n<-length(y)
    if(missing(xreg)) xreg<-cbind(rep(1,n))
    id.ok<-complete.cases(cbind(y,z,xreg))
    y<-y[id.ok]
    z<-z[id.ok]
    xreg<-xreg[id.ok,,drop=FALSE]
    n<-length(y)
    if(family$family=="poisson") size=1
    if(length(size)==1) size<-rep(size,n)
    yN<-y/size
    k=ncol(xreg) #per un modello ~1+x
    if(is.null(values)) values<-seq(min(z), max(z), length=10)
    n1<-length(values)
    PSI<-matrix(values, nrow=n, ncol=n1, byrow=TRUE) #(era X2) matrice di valori di psi
    X1<-matrix(z, nrow=n, ncol=n1, byrow=FALSE) #matrice della variabile Z
    X<-eval(parse(text=fn), list(x=X1, p=PSI)) #X<-pmax(X1-X2,0)
    pmaxMedio<-rowMeans(X)
    
    o<-glm.fit(yN, x=xreg, weights=size, family=family)
    r<-y-(o$fitted*size)
    sc<-drop(crossprod(r, pmaxMedio))
    #    if(output=="unst.score") return(drop(sc))
    
    p <- o$rank
    Qr <- o$qr
    COV <- chol2inv(Qr$qr[1:p, 1:p, drop = FALSE])  #vcov(glm(y~x, family=poisson))
    A<-xreg%*%COV%*%crossprod(xreg, diag(o$weights))
    h<- drop(tcrossprod(pmaxMedio, diag(n)- A))
    if(missing(weights.var)) weights.var<-o$weights
    v.s<- drop(crossprod(h*sqrt(weights.var))) #t(h)%*%diag(exp(lp))%*%h
    
    ris<-if(length(fn)<=1) sc/sqrt(v.s) else drop(crossprod(sc,solve(v.s,sc)))
    #    if(output=="score") return(drop(ris))
    
    pvalue<-  switch(alternative,
                     less = pnorm(ris, lower.tail =TRUE) ,
                     greater = pnorm(ris, lower.tail =FALSE) ,
                     two.sided = 2*pnorm(abs(ris), lower.tail =FALSE) 
    )
    
    #    pvalue<- if(length(fn)<=1) 2*pnorm(abs(ris), lower.tail =FALSE) else pchisq(ris,df=length(fn), lower.tail =FALSE)
    # NB: se calcoli ris<-drop(t(sc)%*%solve(v.s,sc))/(length(fn)*sigma^2) devi usare pf(ris,df1=length(fn),df2=df.t, lower.tail =FALSE)
    return(c(ris, pvalue))
  }
  #----------------------------------------------------
  if(!inherits(obj, "lm")) stop("A '(g)lm', or 'segmented-(g)lm' model is requested")
  fn="pmax(x-p,0)"
  ties.ok=FALSE
  if(missing(seg.Z)){
    if(inherits(obj, "segmented") && length(obj$nameUV$Z)==1) seg.Z<- as.formula(paste("~", obj$nameUV$Z ))
    if(!inherits(obj, "segmented") && length(all.vars(formula(obj)))==2) seg.Z<- as.formula(paste("~", all.vars(formula(obj))[2]))
  } else {
    if(class(seg.Z)!="formula") stop("'seg.Z' should be an one-sided formula")
  }
  if(any(c("$","[") %in% all.names(seg.Z))) stop(" '$' or '[' not allowed in 'seg.Z' ")
  name.Z <- all.vars(seg.Z)
  if(length(name.Z)>1) stop("Only one variable can be specified in 'seg.Z' ")
  
  nomiU.term<-grep(name.Z, obj$nameUV$U, value=TRUE) #termini U per relativi alla variabile nomeZ
  #se length(nomiU.term)==0 la variabile in seg.Z non e' nel modello (si sta assumendo che la left slope ==0)
  if(length(nomiU.term)==0 && more.break) warning(paste("variable", name.Z, "has no breakpoint.. 'more.break=TRUE' ignored"), call.=FALSE)
  #browser()  
  if(k<=1) stop("k>1 requested! k>=10 is recommended")
  if(k<10) warnings("k>=10 is recommended")
  alternative <- match.arg(alternative)
  if(!n.break%in%1:2) stop(" 'n.break' should be 1 or 2", call. = FALSE)
  if(n.break==2) alternative<-"two.sided"
  isGLM<-"glm"%in%class(obj)
  
  if(isGLM){
    if (is.null(dispersion)) dispersion <- summary.glm(obj)$dispersion
    if(inherits(obj, "segmented")){
        if(more.break && !name.Z %in% obj$nameUV$Z) stop(" 'more.break' is meaningful only if at least 1 breakpoint has been estimated")
        Call<-mf<-obj$orig.call #del GLM
        formulaSeg <-formula(obj) #contiene le variabili U e le psi
        formulaNull<- update.formula(formulaSeg, paste("~.-",paste(obj$nameUV$V, collapse="-"))) #rimuovi le variabili "psi.."
        #se length(nomiU.term)==0 la variabile in seg.Z non e' nel modello (si sta assumendo che la left slope ==0)
        if(!more.break && length(nomiU.term)>0){
          if(length(nomiU.term)>1) stop(" 'more.break=FALSE' does not work with multiple breakpoints referring to the same variable specified in seg.Z", call. = FALSE)  
          formulaNull <-update.formula(formulaNull,paste("~.-",paste(nomiU.term, collapse="-")))  
          #non contiene U del termine di interesse, MA contiene eventuali altri termini U
        }
        mf$formula<-formulaNull
        mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+")) #se il modello inziale non contiene seg.Z.. 
        if(!is.null(obj$orig.call$offset) || !is.null(obj$orig.call$weights) || !is.null(obj$orig.call$subset)){ 
          mf$formula <- update.formula(mf$formula, 
                                       paste(".~.+", paste(c(all.vars(obj$orig.call$offset), 
                                                             all.vars(obj$orig.call$weights),
                                                             all.vars(obj$orig.call$subset)), 
                                                           collapse = "+")))
        }
        m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1L]] <- as.name("model.frame")
        for(i in 1:length(obj$nameUV$U)) assign(obj$nameUV$U[i], obj$model[,obj$nameUV$U[i]], envir=parent.frame())
        mf <- eval(mf, parent.frame())
        mt <- attr(mf, "terms")
        #interc<-attr(mt,"intercept")
        y <- model.response(mf, "any")
        X0<- if (!is.empty.model(mt)) model.matrix(mt, mf)
        Z<-X0[ ,match(name.Z, colnames(X0))]
        n<-length(Z)
        if(is.null(values)) values<-seq(min(Z), max(Z), length=k) #values<-seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)
        n1<-length(values)
        X1<-matrix(Z, nrow=n, ncol=n1, byrow=FALSE) #matrice della variabile Z
        PSI<-matrix(values, nrow=n, ncol=n1, byrow=TRUE)
        X<-eval(parse(text=fn), list(x=X1, p=PSI))    #   fn t.c. length(fn)<=1; fn="pmax(x-p,0)" definita sopra..
        pmaxMedio <-as.matrix(rowMeans(X))
        if(n.break>1){
          XX<-sapply(1:length(values), function(.x) X[,-(1:.x), drop=FALSE])
          XX<-do.call("cbind", XX)
          if(ties.ok) XX<-cbind(X, XX)
          pmaxMedio2 <- rowMeans(XX)
          pmaxMedio <- cbind(pmaxMedio, pmaxMedio2)
        }
        if(only.term) return(pmaxMedio)
        #necessario salvare pmaxMedio in mf???
        mf$pmaxMedio<-pmaxMedio 
        Call$formula<- formulaNull
        Call$data<-quote(mf)
        obj0<-eval(Call)
        # pos<-1
        # assign("mf", mf, envir=as.environment(pos))        
        # r<-as.numeric(as.matrix(add1(obj0, ~.+pmaxMedio,  scale=dispersion, test="Rao"))[2,c("scaled Rao sc.", "Pr(>Chi)")])
        ws <- sqrt(obj0$weights[obj0$weights>0])
        res<-obj0$residuals[obj0$weights>0]
        zw <- ws * res 
        A <- qr.resid(obj0$qr, ws * pmaxMedio[obj0$weights>0,])
        u<-t(A)%*% zw
        v<-crossprod(A)
        r<-if(n.break==1) u/sqrt(v*dispersion) else t(u)%*% solve(v) %*%u/dispersion 
        #r<- (colSums(as.matrix(A * zw))/sqrt(colSums(as.matrix(A * A)))/sqrt(dispersion))      
        p2<- if(n.break==1) 2*pnorm(abs(r), lower.tail=FALSE) else pchisq(r, df=n.break, lower.tail=FALSE)
        pvalue<-  switch(alternative,
                         less = pnorm(r, lower.tail =TRUE) ,
                         greater = pnorm(r, lower.tail =FALSE) ,
                         two.sided = p2)
        r<-c(r, pvalue)
        # ================fine se e' GLM+segmented.
    } else { 
      #=================Se e' GLM (non segmented)
      Call<-mf<-obj$call
      mf$formula<-formula(obj)
      m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
      mf <- mf[c(1, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1L]] <- as.name("model.frame")
      formulaNull <- formula(obj)
      mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
      #aggiunto 12/03/18 (non trovava la variable weights perche' era salvata in mf come "(weights)")     
      if(!is.null(obj$call$offset) || !is.null(obj$call$weights) || !is.null(obj$call$subset)){ 
        mf$formula <-update.formula(mf$formula, 
                                    paste(".~.+", paste(
                                      c(all.vars(obj$call$offset), 
                                        all.vars(obj$call$weights),
                                        all.vars(obj$call$subset)), 
                                      collapse = "+")))
      }
      mf <- eval(mf, parent.frame())
      mt <- attr(mf, "terms")
      XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf)
      n <- nrow(XREG)
      Z<- XREG[,match(name.Z, colnames(XREG))]
      if(!name.Z %in% names(coef(obj))) XREG<-XREG[,-match(name.Z, colnames(XREG)),drop=FALSE]
      if(is.null(values)) values<-seq(min(Z), max(Z), length=k) #values<-seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)
      n1<-length(values)
      PSI<-matrix(values, nrow=n, ncol=n1, byrow=TRUE) #(era X2) matrice di valori di psi
      X1<-matrix(Z, nrow=n, ncol=n1, byrow=FALSE) #matrice della variabile Z
      X<-eval(parse(text=fn), list(x=X1, p=PSI))    #   fn t.c.  length(fn)<=1
      pmaxMedio<-as.matrix(rowMeans(X))
      if(n.break>1){
        XX<-sapply(1:length(values), function(.x) X[,-(1:.x), drop=FALSE])
        XX<-do.call("cbind", XX)
        if(ties.ok) XX<-cbind(X, XX)
        pmaxMedio2 <- rowMeans(XX)
        pmaxMedio <- cbind(pmaxMedio, pmaxMedio2)
      }
      if(only.term) return(pmaxMedio)
      #r<-as.numeric(as.matrix(add1(update(obj, data=mf), ~.+pmaxMedio,  scale=dispersion, test="Rao"))[2,c("scaled Rao sc.", "Pr(>Chi)")])       
      #Call$formula<- formulaNull
      #Call$data<-quote(mf)
      #obj0<-eval(Call)
      ws <- sqrt(obj$weights[obj$weights>0])
      res<-obj$residuals[obj$weights>0]
      zw <- ws * res
      A <- qr.resid(obj$qr, ws * pmaxMedio[obj$weights>0,])
      u<-t(A)%*% zw
      v<-crossprod(A)
      r<-if(n.break==1) u/sqrt(v*dispersion) else t(u)%*% solve(v) %*%u/dispersion 
      #r<- (colSums(as.matrix(A * zw))/sqrt(colSums(as.matrix(A * A)))/sqrt(dispersion))      
      p2<- if(n.break==1) 2*pnorm(abs(r), lower.tail=FALSE) else pchisq(r, df=n.break, lower.tail=FALSE)
      pvalue<-  switch(alternative,
                       less = pnorm(r, lower.tail =TRUE) ,
                       greater = pnorm(r, lower.tail =FALSE) ,
                       two.sided = p2)
      r<-c(r, pvalue)
      
      #fine se e' un GLM
    }
    
  } else { ##============================== Se e' un LM..
    if(is.null(dispersion)) dispersion<- summary(obj)$sigma^2
    if(is.null(df.t)) df.t <- obj$df.residual
    #df.ok<- if(!is.null(df.t)) df.t else obj$df.residual
    
    #============ se e' LM+segmented
    if(inherits(obj, "segmented")){ 
      if(more.break && !name.Z %in% obj$nameUV$Z) stop(" stop 'more.break' is meaningful only if at least 1 breakpoint has been estimated", call.=FALSE )
      Call<-mf<-obj$orig.call
      mf$formula<-formula(obj)
      m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
      mf <- mf[c(1, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1L]] <- as.name("model.frame")
      mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
      #formulaOrig<-formula(obj)
      if(!is.null(obj$orig.call$offset) || !is.null(obj$orig.call$weights) || !is.null(obj$orig.call$subset)){ 
        mf$formula <- update.formula(mf$formula, 
                                     paste(".~.+", paste(c(all.vars(obj$orig.call$offset), 
                                                           all.vars(obj$orig.call$weights),
                                                           all.vars(obj$orig.call$subset)), 
                                                         collapse = "+")))
      }
      mf$formula<-update.formula(mf$formula,paste("~.-",paste(obj$nameUV$V, collapse="-"))) #rimuovi le variabili "psi.."
      if(!more.break) {
        if(length(nomiU.term)>1) stop(" 'more.break=FALSE' does not work with multiple breakpoints referring to the same variable specified in seg.Z", call. = FALSE)
        #ovvero il test funziona per un solo breakpoint..
        mf$formula<-update.formula(mf$formula,paste("~.-",paste(nomiU.term, collapse="-"))) #rimuovi il termine U in questione, cioe' solo per una variabile
        #altre variabili "U" relative a piu' variabili devono rimanere.. 
      }
      formulaNull <- formula(mf)
      
      ###############
      #PERCHE' NON estrarre direttamente la model.matrix(obj)
      #X <-model.matrix(obj)
      #X <- X[, !(colnames(X) %in% obj$nameUV$V), drop=FALSE]
      #perche' poi se la variabile NON e' nel modello (perche' left slope=0) non so come recuperarla
      
      if(more.break) {
        mf$data<-quote(model.frame(obj))
        mf<-eval(mf)
        } else {
        mf <- eval(mf, parent.frame())  
      }
      #for(i in 1:length(obj$nameUV$U)) assign(obj$nameUV$U[i], obj$model[,obj$nameUV$U[i]], envir=parent.frame())
      
      y <- model.response(mf, "any")
      weights <- as.vector(model.weights(mf))
      offset <- as.vector(model.offset(mf))
      
      mt <- attr(mf, "terms")
      #interc<-attr(mt,"intercept")
      
      X0<- if (!is.empty.model(mt)) model.matrix(mt, mf)
      Z<-X0[ ,match(name.Z, colnames(X0))]
      n<-length(Z)
      if(is.null(values)) values<-seq(min(Z), max(Z), length=k) #values<-seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)
#      n1<-length(values)
#      X1<-matrix(Z, nrow=n, ncol=n1, byrow=FALSE) #matrice della variabile Z
#      PSI<-matrix(values, nrow=n, ncol=n1, byrow=TRUE)
#      X<-eval(parse(text=fn), list(x=X1, p=PSI))    #   fn t.c. length(fn)<=1; fn="pmax(x-p,0)" definita sopra..
#      mf$pmaxMedio<- pmaxMedio <-rowMeans(X)
      
      r<-test.Sc2(y=y, z=Z, xreg=X0, sigma=sqrt(dispersion), values=values, fn=fn, df.t=df.t, alternative=alternative, 
                  w=weights, offs=offset, nbreaks=n.break, ties.ok=FALSE, only.term=only.term)
      
      #fine se e' LM+segmented. 
    } else {
      #=================Se e' LM (non segmented)
      Call<-mf<-obj$call
      mf$formula<-formula(obj)
      m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
      mf <- mf[c(1, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1L]] <- as.name("model.frame")
      formulaNull <- formula(obj)
      mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
      #aggiunto 12/03/18 (non trovava la variable weights perche' era salvata in mf come "(weights)")     
      if(!is.null(obj$call$offset) || !is.null(obj$call$weights) || !is.null(obj$call$subset)){ 
        mf$formula <-update.formula(mf$formula, 
                                    paste(".~.+", paste(
                                      c(all.vars(obj$call$offset), 
                                        all.vars(obj$call$weights),
                                        all.vars(obj$call$subset)), 
                                      collapse = "+")))
      }
      mf <- eval(mf, parent.frame())
      y <- model.response(mf, "any")
      weights <- as.vector(model.weights(mf))
      offset <- as.vector(model.offset(mf))
      
      mt <- attr(mf, "terms")
      XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf)
      n <- nrow(XREG)
      Z<- XREG[,match(name.Z, colnames(XREG))]
      if(!name.Z %in% names(coef(obj))) XREG<-XREG[,-match(name.Z, colnames(XREG)),drop=FALSE]
      if(is.null(values)) values<-seq(min(Z), max(Z), length=k) #values<-seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)
      #n1<-length(values)
      #PSI<-matrix(values, nrow=n, ncol=n1, byrow=TRUE) #(era X2) matrice di valori di psi
      #X1<-matrix(Z, nrow=n, ncol=n1, byrow=FALSE) #matrice della variabile Z
      #X<-eval(parse(text=fn), list(x=X1, p=PSI))    #   fn t.c.  length(fn)<=1
      #pmaxMedio<-rowMeans(X)
      r<-test.Sc2(y=y, z=Z, xreg=XREG, sigma=sqrt(dispersion), values=values, fn=fn, df.t=df.t, alternative=alternative, 
                  w=weights, offs=offset, nbreaks=n.break, ties.ok=FALSE, only.term=only.term)
      #r<-as.numeric(as.matrix(add1(update(obj, data=mf), ~.+pmaxMedio,  scale=dispersion, test="Rao"))[2,c("scaled Rao sc.", "Pr(>Chi)")])       
    }
  } #end se LM
  #################################################
  if(only.term) return(r)
  #################################################
  if(is.null(obj$family$family)) {
    famiglia<-"gaussian"
    legame<-"identity"
  } else {
    famiglia<-obj$family$family
    legame<-obj$family$link
  }
  
  out <- list(method = "Score test for one/two changes in the slope",
              data.name=paste("formula =", as.expression(formulaNull), "\nbreakpoint for variable =", name.Z, 
                              "\nmodel =",famiglia,", link =", legame ,", method =", obj$call[[1]]),
              statistic = c(`observed value` = r[1]),
              parameter = c(n.points = length(values)), p.value = r[2],
              #alternative = paste(alternative, " (",n.break ,"breakpoint ) ")
              alternative = paste(alternative,"   (",n.break ,if(n.break==1) " breakpoint) " else " breakpoints) ", sep=""))
  class(out) <- "htest"
  return(out)
}



