`davies.test` <- 
function (obj, seg.Z, k = 10, alternative = c("two.sided", "less", "greater"), beta0=0, dispersion=NULL) {
    extract.t.value.U<-function(x){
        #estrae il t-value dell'ultimo coeff in un oggetto restituito da lm.fit
            #x<-x$obj
            R<-qr.R(x$qr)
            p<-ncol(R)
            n<-length(x$fitte.values)
            invR<-backsolve(R,diag(p))
            hat.sigma2<-sum(x$residuals^2)/(n-p)
            #solve(crossprod(qr.X(x$qr)))
            V<-tcrossprod(invR)*hat.sigma2
            tt<-x$coefficients[p]/sqrt(V[p,p])
            tt}
    extract.t.value.U.glm<-function(object,dispersion,isGLM,beta0){
        #estrae il t-value dell'ultimo coeff in un oggetto restituito da lm.wfit/glm.fit
        est.disp <- FALSE
        df.r <- object$df.residual
        if (is.null(dispersion))
          dispersion <- if(isGLM&&(object$family$family%in%c("poisson","binomial"))) 1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0))
                warning("observations with zero weight not used for calculating dispersion")
            sum((object$weights * object$residuals^2)[object$weights > 0])/df.r
        }  else {
            est.disp <- TRUE
            NaN
        }
    p <- object$rank
    p1 <- 1L:p
    Qr <- object$qr
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
    covmat <- dispersion * covmat.unscaled
    tvalue <- (coef.p[p]-beta0)/sqrt(covmat[p,p])
    tvalue
    }#end extract.t.value.U.glm 

    
    alternative <- match.arg(alternative)
#-------------------------------------------------------------------------------
#    it.max <- old.it.max<- 0
#    toll <- .1
#    visual <- FALSE
#    stop.if.error<-TRUE
#    last <- FALSE
#    K<-10
#    h<-min(abs(control$h),1)
#
    if(length(all.vars(seg.Z))>1) warning("multiple segmented variables ignored in 'seg.Z'",call.=FALSE)
    isGLM<-"glm"%in%class(obj)
    Call<-mf<-obj$call
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
    mf <- eval(mf, parent.frame())
    weights <- as.vector(model.weights(mf))
    offs <- as.vector(model.offset(mf))
    if(!is.null(Call$weights)){ #"(weights)"%in%names(mf)
      names(mf)[which(names(mf)=="(weights)")]<-as.character(Call$weights)
      #aggiungere???
      # mf["(weights)"]<-weights
      }
    mt <- attr(mf, "terms")
    interc<-attr(mt,"intercept")
    y <- model.response(mf, "any")
    XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    n <- nrow(XREG)
#    weights <- as.vector(model.weights(mf))
#    offs <- as.vector(model.offset(mf))
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(offs)) offs <- rep(0, n)

    id.duplic<-match(all.vars(formula(obj)),all.vars(seg.Z),nomatch=0)>0
    if(any(id.duplic)) {
        #new.mf<-mf[,id.duplic,drop=FALSE]
        new.mf<-mf[,all.vars(formula(obj))[id.duplic],drop=FALSE]
#        new.XREGseg<-model.matrix(attr(new.mf, "terms"), new.mf, contrasts)
        new.XREGseg<-data.matrix(new.mf)
        XREG<-cbind(XREG,new.XREGseg)
        }
    n.Seg <- 1
    #n.Seg <- if(is.list(psi)) length(psi) else 1
    #n.psi<- length(unlist(psi))
    n.psi<- 1
    id.n.Seg<-(ncol(XREG)-n.Seg+1):ncol(XREG)
    XREGseg<-XREG[,id.n.Seg,drop=FALSE]
    XREG<-XREG[,match(c("(Intercept)",all.vars(formula(obj))[-1]),colnames(XREG),nomatch =0),drop=FALSE]
    
    
    n <- nrow(XREG)
    Z<-lapply(apply(XREGseg,2,list),unlist) #prende anche i nomi!
    name.Z <- names(Z) <- colnames(XREGseg)
#    if(length(Z)>1) stop("Only one single segmented variable is allowed")

    dev0<-if(isGLM) obj$dev else sum(obj$residuals^2)
    nomiU <- paste("U1", name.Z, sep = ".")
    nomiOK<-nomiU
    opz<-list(toll=.1,h=1,stop.if.error=FALSE,dev0=dev0,visual=FALSE,it.max=0,nomiOK=nomiOK)
    Z<-matrix(unlist(Z),ncol=1)
    qq <- quantile(Z, prob = c(0.05, 0.95), names = FALSE, na.rm = TRUE)
    sx <- qq[1]
    dx <- qq[2]
    valori <- seq(sx, dx, length = k)
    ris.valori <- NULL
    eta0<-obj$linear.predictors
    
    for (i in 1:length(valori)) {
        psi<-valori[i]
        XX<-cbind(XREG, pmax((Z - psi), 0))
        #PSI <- matrix(rep(psi, n), ncol = 1)
        #ob <- segmented.lm.fit(y,XREG,Z,PSI,w,o,opz)
        ob<-if(!isGLM) lm.wfit(x = XX, y = y, w = weights, offset = offs)
              else glm.fit(x = XX, y = y, weights = weights, offset = offs, 
                  family=family(obj), etastart=eta0)        
        eta0<-ob$linear.predictors
        ris.valori[(length(ris.valori)) + 1] <- if (is.list(ob)) {
            #summary(o)$coef[length(coef(ogg)) + 1, 3]
        #extract.t.value.U(ob)
        extract.t.value.U.glm(ob,dispersion,isGLM,beta0)
        } else { NA }
    }
    valori <- valori[!is.na(ris.valori)]
    ris.valori <- ris.valori[!is.na(ris.valori)]
    V <- sum(abs(diff(ris.valori)))

    onesided <- TRUE
    if (alternative == "less") {
        M <- min(ris.valori)
        best<-valori[which.min(ris.valori)]
        p.naiv <- pnorm(M, lower = TRUE)
    }
    else if (alternative == "greater") {
        M <- max(ris.valori)
        best<-valori[which.max(ris.valori)]
        p.naiv <- pnorm(M, lower = FALSE)
    }
    else {
        M <- max(abs(ris.valori))
        best<-valori[which.max(abs(ris.valori))]
        p.naiv <- pnorm(M, lower = FALSE)
        onesided <- FALSE
    }
    approxx <- V * exp(-(M^2)/2)/sqrt(8 * pi)
    p.adj <- p.naiv + approxx
    p.adj <- ifelse(onesided, 1, 2) * p.adj
    if(is.null(obj$family$family)) {
          famiglia<-"gaussian"
          legame<-"identity"} else {
               famiglia<-obj$family$family
               legame<-obj$family$link
          }
    out <- list(method = "Davies' test for a change in the slope",
        data.name=paste("Model = ",famiglia,", link =", legame,
        "\nformula =", as.expression(formula(obj)),
#        data.name = paste(as.expression(ogg$call),
        "\nsegmented variable =", name.Z),
        statistic = c("`Best' at" = best),
        parameter = c(n.points = length(valori)), p.value = min(p.adj,1),
        alternative = alternative)
    class(out) <- "htest"
    return(out)
}
    
    
