`pscore.test` <- function(obj, seg.Z, k = 10, alternative = c("two.sided", "less", "greater"),
    values=NULL, dispersion=NULL, df.t=NULL, more.break=FALSE) { 
#-------------------------------------------------------------------------------
test.Sc2<-function(y, z, xreg, sigma=NULL, values=NULL, fn="pmax(x-p,0)", df.t="Inf", alternative, w=NULL){
#xreg: la matrice del disegno del modello nullo. Se mancante viene assunta solo l'intercetta.
#Attenzione che se invXtX e xx vengono entrambe fornite, non viene fatto alcun controllo
#invXtX: {X'X}^{-1}. if missing it is computed from xreg
#sigma: the sd. If missing it is computed from data (under the *null* model)
#values: the values with respect to ones to compute the average term. If NULL 10 values from min(z) to max(z) are taken.
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
          X<-eval(parse(text=fn), list(x=X1, p=PSI)) #X<-pmax(X1-X2,0)
          pmaxMedio<-rowMeans(X)
          } else {
          pmaxMedio<-matrix(,n,length(fn))
          #list.X<-vector("list", length=length(fn))
          for(j in 1:length(fn)){
              #list.X[[j]]<-eval(parse(text=fn[j]), list(x=X1, p=PSI))
              X<-eval(parse(text=fn[[j]]), list(x=X1, p=PSI))
              pmaxMedio[,j]<-rowMeans(X)
              }
          }
       }
       if(is.null(w)){
              invXtX<-solve(crossprod(xreg))
              IA<-diag(n) - xreg%*%tcrossprod(invXtX, xreg) #I- hat matrix
              r<-IA%*%y
              sc<- sum(pmaxMedio*r)
              v.s<-crossprod(pmaxMedio,IA)%*% pmaxMedio
       } else {
              invXtX<-solve(crossprod(sqrt(w)*xreg))
              IA<-diag(n) - xreg%*%tcrossprod(invXtX, xreg*w) #I-hat matrix
              sc<-t(pmaxMedio*w) %*% IA  %*% y
              v.s<- t(pmaxMedio*w) %*% crossprod(t(IA)/sqrt(w))%*%(w*pmaxMedio)
       }
       
#       ris<-if(length(fn)<=1) sc/(sigma*sqrt(v.s)) else drop(crossprod(sc,solve(v.s,sc)))/(sigma^2)
#       if(length(fn)<=1 && cadj) ris<- sign(ris)*sqrt((ris^2)*(1-(3-(ris^2))/(2*n)))
        ris<- drop(sc/(sigma*sqrt(v.s))) 
       df.t<-eval(parse(text=df.t))
       
       pvalue<-  switch(alternative,
         less = pt(ris, df=df.t, lower.tail =TRUE) ,
         greater = pt(abs(ris), df=df.t, lower.tail =FALSE) ,
         two.sided = 2*pt(abs(ris), df=df.t, lower.tail =FALSE) 
         )
       #pvalue<- 2*pt(abs(ris), df=df.t, lower.tail =FALSE) 
       r<-c(ris, pvalue, pmaxMedio)
       r
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
         greater = pnorm(abs(ris), lower.tail =FALSE) ,
         two.sided = 2*pnorm(abs(ris), lower.tail =FALSE) 
         )
    
#    pvalue<- if(length(fn)<=1) 2*pnorm(abs(ris), lower.tail =FALSE) else pchisq(ris,df=length(fn), lower.tail =FALSE)
    # NB: se calcoli ris<-drop(t(sc)%*%solve(v.s,sc))/(length(fn)*sigma^2) devi usare pf(ris,df1=length(fn),df2=df.t, lower.tail =FALSE)
    return(c(ris, pvalue))
    }
#----------------------------------------------------

    fn="pmax(x-p,0)"
    if(inherits(obj, "glm")) stop("Currently only 'lm', or 'segmented-lm' models are allowed")
    if(!inherits(obj, "lm")) stop("A 'lm', or 'segmented-lm' model is requested")
    if(inherits(obj, "segmented") && length(obj$nameUV$Z)==1) seg.Z<- as.formula(paste("~", obj$nameUV$Z ))
    if(!inherits(obj, "segmented") && length(all.vars(formula(obj)))==2) seg.Z<- as.formula(paste("~", all.vars(formula(obj))[2]))
    if(class(seg.Z)!="formula") stop("'seg.Z' should be an one-sided formula")
    name.Z <- all.vars(seg.Z)
    if(length(name.Z)>1) stop("Only a single segmented variable can be specified in 'seg.Z' ")
    
    if(k<=1) stop("k>1 requested! k>=10 is recommended")
    if(k<10) warnings("k>=10 is recommended")
    alternative <- match.arg(alternative)
    #if(length(all.vars(seg.Z))>1) warning("multiple segmented variables ignored in 'seg.Z'",call.=FALSE)
    isGLM<-"glm"%in%class(obj)
    
    if(inherits(obj, "segmented")){
        if(is.null(dispersion)) {
        dispersion<- if(inherits(obj, "glm")) summary(obj)$dispersion else summary(obj)$sigma^2   
        }
        df.ok<- if(!is.null(df.t)) df.t else obj$df.residual
        if(!is.null(eval(obj$call$obj)$call$data)) mf$data <- eval(obj$call$obj)$call$data
        y<- obj$res+obj$fitted        
        if(!is.null(obj$offset)) y<- y-obj$offset
        weights<- obj$weights
        X0<-model.matrix(obj)
        Z<-X0[,name.Z]
        #escludi *tutte* le variabili psi (sia da X0 che dalla formula)
        X0<-X0[, -match(obj$nameUV$V, colnames(X0))]
        formulaNull <-update.formula(formula(obj),paste("~.-",paste(obj$nameUV$V, collapse="-"))) #togli tutti i termini "V"
        #se vuoi saggiare per un additional breakpoint
        if(!more.break){
          idU<-grep(name.Z, obj$nameUV$U)
          X0<-X0[, -match(obj$nameUV$U[idU], colnames(X0))]
          formulaNull <-update.formula(formulaNull, paste("~.-",paste(obj$nameUV$U[idU], collapse="-")))
          }
        #browser()
        if(is.null(values)) values<-seq(min(Z), max(Z), length=k) #values<-seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)
#        formulaOrig<-mf$formula<-update.formula(mf$formula,paste("~.-",paste(obj$nameUV$V, collapse="-"))) #togli tutti i termini "V"
#        update.formula(formulaOrig, paste("~.-",paste(obj$nameUV$U[idU], collapse="-")))
#        #for(i in 1:length(obj$nameUV$U)) assign(obj$nameUV$U[i], obj$model[,obj$nameUV$U[i]], envir=parent.frame())
#        formulaOrig<-update.formula(formulaOrig, paste("~.-",paste(obj$nameUV$U, collapse="-")))
         r<-test.Sc2(y=y, z=Z, xreg=X0, sigma=sqrt(dispersion), values=values, fn=fn, df.t=df.ok, alternative=alternative, w=weights)
        } else {
        #SE l'oggetto obj NON e' "segmented"..
        # ci possono essere mancanti nella variabile seg.Z, quindi devi costruire un nuovo dataframe....    
        Call<-mf<-obj$call
        mf$formula<-formula(obj)
        m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1L]] <- as.name("model.frame")
        mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
        formulaNull <- formula(obj)
        mf <- eval(mf, parent.frame())

        weights <- as.vector(model.weights(mf))
        offs <- as.vector(model.offset(mf))
        
       mt <- attr(mf, "terms")
       interc<-attr(mt,"intercept")
       y <- model.response(mf, "any")
       XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
       n <- nrow(XREG)
       if (is.null(weights)) weights <- rep(1, n)
       if (!is.null(offs)) y<-y-offs 
       Z<- XREG[,match(name.Z, colnames(XREG))]
       if(!name.Z %in% names(coef(obj))) XREG<-XREG[,-match(name.Z, colnames(XREG)),drop=FALSE]
       if(is.null(values)) values<-seq(min(Z), max(Z), length=k) #values<-seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)
       if(is.null(dispersion)) {
         dispersion<- if(class(obj)[1]=="glm") sqrt(summary.glm(obj)$dispersion) else summary(obj)$sigma^2
         df.ok<- if(!is.null(df.t)) df.t else obj$df.residual
       }

       

       ##---------se obj e' modello nullo (non contiene la variabile segmented):
#questo dovrebbe funzionare per anche GLM/LM (attenzione che il p-valore e' calcolato da add1() con la z e non t )
#    n1<-length(values)
#    PSI<-matrix(values, nrow=n, ncol=n1, byrow=TRUE) #(era X2) matrice di valori di psi
#    X1<-matrix(Z, nrow=n, ncol=n1, byrow=FALSE) #matrice della variabile Z
#    X<-eval(parse(text=fn), list(x=X1, p=PSI)) #      fn t.c.  length(fn)<=1
#    pmaxMedio<-rowMeans(X)
#    #mf<-model.frame(obj)
#    mf$pmaxMedio <- pmaxMedio
#    if(!is.null(weights)){ #"(weights)"%in%names(mf)
#            names(mf)[which(names(mf)=="(weights)")]<-all.vars(Call$weights) #as.character(Call$weights)
#            mf["(weights)"]<-weights
#        }
#    r<-as.numeric(add1(update(obj, data=mf), ~.+pmaxMedio,  scale=sigma^2, test="Rao")[c("scaled Rao sc.", "Pr(>Chi)")][2,])
#     r[2]<-1-pf(r[1],df1=1, df2=df.ok)

    r<-test.Sc2(y=y, z=Z, xreg=XREG, sigma=sqrt(dispersion), values=values, fn=fn, df.t=df.ok, alternative=alternative, w=weights)
    }
    
        if(is.null(obj$family$family)) {
          famiglia<-"gaussian"
          legame<-"identity"
            } else {
               famiglia<-obj$family$family
               legame<-obj$family$link
              }

    out <- list(method = "Score test for one change in the slope",
#        data.name=paste("Model = ",famiglia,", link =", legame,
#        "\nformula =", as.expression(formulaOrig),
#        "\nsegmented variable =", name.Z),
        data.name=paste("formula =", as.expression(formulaNull), ",   method =", obj$call[[1]] ,
        "\nmodel =",famiglia,", link =", legame, NULL ,
        "\nsegmented variable =", name.Z),
        statistic = c(`observed value` = r[1]),
        parameter = c(n.points = length(values)), p.value = r[2],
        alternative = alternative)
    class(out) <- "htest"
    return(out)
}








