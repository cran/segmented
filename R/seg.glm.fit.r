seg.glm.fit<-function(y,XREG,Z,PSI,w,o,opz){
#opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,dev0=dev0)
    eta0<-opz$eta0
    fam<-opz$fam
    maxit.glm<-opz$maxit.glm
    #--------------
    nomiOK<-opz$nomiOK
    toll<-opz$toll
    h<-opz$h
    stop.if.error<-opz$stop.if.error
    dev.new<-opz$dev0
    visual<-opz$visual
    it.max<-old.it.max<-opz$it.max
    rangeZ <- apply(Z, 2, range)
    #k<-ncol(Z)
    psi<-PSI[1,]
    H<-1
    it <- 1
    epsilon <- 10
    psi.values <- NULL
    while (abs(epsilon) > toll) {
        k<-ncol(Z)
        U <- pmax((Z - PSI), 0)
        V <- ifelse((Z > PSI), -1, 0)
        #dev.old <- sum(obj$residuals^2)
        X <- cbind(XREG, U, V)
        rownames(X) <- NULL
        if (ncol(V) == 1)
            colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c("U", "V")
        else colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U",
            1:ncol(U), sep = ""), paste("V", 1:ncol(V), sep = ""))
        #obj <- lm.wfit(x = X, y = y, w = w, offset = o)
        #controlla******************
        obj <- suppressWarnings(glm.fit(x = X, y = y, offset = o,
            weights = w, family = fam, control = glm.control(maxit = maxit.glm),
            etastart = eta0))
        eta0 <- obj$linear.predictors
        #---
        dev.old<-dev.new
        dev.new <- obj$dev
        if (visual) {
            if (it == 1)
                cat(0, " ", formatC(dev.old, 3, format = "f"),
                  "", "(No breakpoint(s))", "\n")
            spp <- if (it < 10) "" else NULL
            cat(it, spp, "", formatC(dev.new, 3, format = "f"),
                "\n")
        }
        epsilon <- (dev.new - dev.old)/dev.old
        obj$epsilon <- epsilon
        it <- it + 1
        obj$it <- it
        #class(obj) <- c("segmented", class(obj))
        #list.obj[[length(list.obj) + ifelse(last == TRUE, 0, 1)]] <- obj
        if (k == 1) {
            beta.c <- coef(obj)["U"]
            gamma.c <- coef(obj)["V"]
        }
        else {
            beta.c <- coef(obj)[paste("U", 1:ncol(U), sep = "")]
            gamma.c <- coef(obj)[paste("V", 1:ncol(V), sep = "")]
        }
        if (it > it.max) break
        psi.values[[length(psi.values) + 1]] <- psi.old <- psi
        if(it>=old.it.max && h<1) H<-h
        psi <- psi.old + H*gamma.c/beta.c
        PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
        #check if psi is admissible..
        a <- apply((Z <= PSI), 2, all) #prima era solo <
        b <- apply((Z >= PSI), 2, all) #prima era solo >
        if(stop.if.error) {
            if(sum(a + b) != 0 || is.na(sum(a + b))) stop("(Some) estimated psi out of its range")
            } else {
            id.psi.ok<-!is.na((a+b)<=0)&(a+b)<=0
            Z <- Z[,id.psi.ok,drop=FALSE]
            psi <- psi[id.psi.ok]
            PSI <- PSI[,id.psi.ok,drop=FALSE]
            nomiOK<-nomiOK[id.psi.ok] #salva i nomi delle U per i psi ammissibili
            if(ncol(PSI)<=0) return(0)
            } #end else
        obj$psi <- psi
    } #end while
    obj<-list(obj=obj,it=it,psi=psi,psi.values=psi.values,U=U,V=V,rangeZ=rangeZ,
        epsilon=epsilon,nomiOK=nomiOK)
    return(obj)
    }
