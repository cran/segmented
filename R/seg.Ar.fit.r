seg.Ar.fit<-function(obj, XREG, Z, PSI, opz, return.all.sol=FALSE){
#-----------------
    useExp.k=TRUE
    #-----------------
    est.k<-function(x1,y1,L0){
        ax<-log(x1)
        .x<-cbind(1,ax,ax^2)
        b<-drop(solve(crossprod(.x),crossprod(.x,y1)))
        const<-b[1]-L0
        DD<-sqrt(b[2]^2-4*const*b[3])
        kk<-exp((-b[2]+ DD) /(2*b[3]))
        return(round(kk))
        
        #  ff<-function(xx) b[1]+b[2]*xx + b[3]*xx^2+ L0
        #  a<-uniroot(ff, c(log(x[4]), 3.4))
    }
    #-----------------
    dpmax<-function(x,y,pow=1){
        #deriv pmax
        if(pow==1) -(x>y) #ifelse(x>y, -1, 0)
        else -pow*((x-y)*(x>y))^(pow-1)#-pow*pmax(x-y,0)^(pow-1)
    }
    #-----------
    in.psi<-function(LIM, PSI, ret.id=TRUE){
        #check if psi is inside the range
        a<-PSI[1,]<=LIM[1,]
        b<-PSI[1,]>=LIM[2,]
        is.ok<- !a & !b #TRUE se psi e' OK
        if(ret.id) return(is.ok)
        isOK<- all(is.ok) && all(!is.na(is.ok))
        isOK}
    #------------
    far.psi<-function(Z, PSI, id.psi.group, ret.id=TRUE, fc=.93) {
        #check if psi are far from the boundaries ..s
        #   returns TRUE, if fine.
        #id.far.ok<-sapply(unique(id.psi.group), function(.x) (table(rowSums(((Z>PSI)[,id.psi.group==.x,drop=FALSE])))>=2)[-1]) #[-1] esclude lo zero, x<psi[1] 
        #id.far.ok<-sapply(unique(id.psi.group), function(.x) (tabulate(rowSums(((Z>PSI)[,id.psi.group==.x,drop=FALSE]))+1)>=2)[-1]) #[-1] esclude lo zero, x<psi[1]
        #16/01/20:
        # se un psi assume l'estremo superiore "Z>PSI" non se ne accorge, mentre Z>=PSI, si.. Il contrario e vero con estremo inf e Z>PSI
        nSeg<-length(unique(id.psi.group))
        npsij<-tapply(id.psi.group,id.psi.group,length)
        nj<-sapply(unique(id.psi.group), function(.x) { tabulate(rowSums((Z>PSI)[,id.psi.group==.x,drop=FALSE])+1) }, simplify = FALSE)    
        ff<-id.far.ok<-vector("list",length=nSeg) 
        for(i in 1:nSeg){
            if(length(nj[[i]])!=npsij[i]+1) nj[[i]]<- tabulate(rowSums((Z>=PSI)[,id.psi.group==i,drop=FALSE])+1)
            id.ok<-(nj[[i]] >= 2)
            id.far.ok[[i]] <- id.ok[-length(id.ok)] & id.ok[-1] #esattamente uguale al numero di psi del gruppo i
            ff[[i]]<-ifelse(diff(nj[[i]])>0, 1/fc, fc)
        }
        id.far.ok<-unlist(id.far.ok)
        ff<-unlist(ff)
        if(!ret.id) {return(all(id.far.ok))
        } else {
            attr(id.far.ok,"factor") <- ff
            return(id.far.ok) 
        }
        #if(ret.id) return(id.far.ok) else return(all(id.far.ok))
    } #end far.psi
    #-----------
    adj.psi<-function(psii, LIM) {pmin(pmax(LIM[1,],psii),LIM[2,])} 
    #-----------
    
    n<-nrow(Z)#length(y)
    min.step<-opz$min.step
    rangeZ <- apply(Z, 2, range)
    alpha<-opz$alpha
    limZ <- apply(Z, 2, quantile, names=FALSE, probs=c(alpha,1-alpha))
    
    digits<-opz$digits
    pow<-opz$pow
    nomiOK<-opz$nomiOK
    toll<-opz$toll
    h<-opz$h
    conv.psi<-opz$conv.psi 
    gap<-opz$gap
    #stop.if.error<-opz$stop.if.error
    #fix.npsi<-opz$fix.npsi
    fix.npsi<-opz$stop.if.error
    dev.new<-opz$dev0
    visual<-opz$visual
    id.psi.group<-opz$id.psi.group
    it.max<-old.it.max<-opz$it.max
    psi<-PSI[1,]
    fc<-opz$fc
    names(psi)<-id.psi.group
    epsilon <- 10
    dev.values<-psi.values <- NULL
    #id.psi.ok<-rep(TRUE, length(psi))
    
    it <- 0
    epsilon <- 10
    k.values<-dev.values<- NULL
    psi.values <-list()
    psi.values[[length(psi.values) + 1]] <- NA
    
    
    nomiU<- opz$nomiU
    nomiV<- opz$nomiV
    call.ok <- opz$call.ok   #la call con U e V
    call.ok.noinit<-call.ok
    call.ok.noinit$init<-NULL
    
    call.noV <- opz$call.noV #la call con U ONLY
    call.noV.noinit<-call.noV
    call.noV.noinit$init<-NULL
    
    toll<-opz$toll
    k<-ncol(Z)
    mio.init<-NULL
    mio.init.noV<-NULL
    if(!in.psi(limZ,PSI,FALSE))  stop("starting psi out of the range", call.=FALSE)
    if(!far.psi(Z,PSI,id.psi.group,FALSE)) stop("psi values too close each other. Please change (decreases number of) starting values", call.=FALSE)
    n.psi1<-ncol(Z)
    #==============================================
    U <- ((Z-PSI)*(Z>PSI)) #pmax((Z - PSI), 0)^pow[1]
    colnames(U)<-nomiU
    if(pow[1]!=1) U<-U^pow[1]
    obj0 <- suppressWarnings(try(eval(call.noV.noinit), silent=TRUE))
    if(class(obj0)[1]=="try-error")  stop("The first fit with U variables does not work..", call.=FALSE)#obj0 <- suppressWarnings(eval(call.noV.noinit)) ##a volte con i valori iniziali arima() )non converge!! Quindi provo senza init
    L0<- -obj0$loglik
    n.intDev0<-nchar(strsplit(as.character(L0),"\\.")[[1]][1])
    dev.values[length(dev.values) + 1] <- opz$dev0 #del modello iniziale (senza psi)
    dev.values[length(dev.values) + 1] <- L0 #modello con psi iniziali
    psi.values[[length(psi.values) + 1]] <- psi #psi iniziali
    #==============================================
    if (visual) {
        cat(paste("iter = ", sprintf("%2.0f",0),
                  "  llik = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), -L0), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
                  "  k = ", sprintf("%2.0f", NA),
                  "  n.psi = ",formatC(length(unlist(psi)),digits=0,format="f"), 
                  "  ini.psi = ",paste(formatC(unlist(psi),digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                  sep=""), "\n")
    }
    #==============================================  
    id.warn <- FALSE
    id.psi.changed<-rep(FALSE, it.max)
    while (abs(epsilon) > toll) {
        it<-it+1
        n.psi0 <- n.psi1
        n.psi1 <- ncol(Z)
        if(n.psi1!=n.psi0){
            U <- ((Z-PSI)*(Z>PSI)) #pmax((Z - PSI), 0)^pow[1]
            if(pow[1]!=1) U<-U^pow[1]
            obj0 <- suppressWarnings(try(eval(call.noV), silent=TRUE))
            if(class(obj0)[1]=="try-error")  obj0 <- suppressWarnings(eval(call.noV.noinit)) ##a volte con i valori iniziali arima() )non converge!! Quindi provo senza init
            L0<- -obj0$loglik
        } 
        V <- dpmax(Z,PSI,pow=pow[2])# ifelse((Z > PSI), -1, 0)
        X <- cbind(XREG, U, V)
        rownames(X) <- NULL
        colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U", 1:ncol(U), sep = ""), paste("V", 1:ncol(V), sep = ""))
        obj <- suppressWarnings(try(eval(call.ok), silent=TRUE))
        if(class(obj)[1]=="try-error")  obj <- suppressWarnings(eval(call.ok.noinit)) ##a volte con i valori iniziali arima() )non converge!! Quindi provo senza init
        beta.c <- coef(obj)[paste("U", 1:ncol(U), sep = "")]
        gamma.c <- coef(obj)[paste("V", 1:ncol(V), sep = "")]

        if(any(is.na(c(beta.c, gamma.c)))){
            if(fix.npsi)  {
                if(return.all.sol) return(list(dev.values, psi.values)) else stop("breakpoint estimate too close or at the boundary causing NA estimates.. too many breakpoints being estimated?", call.=FALSE)
            } else {
                id.coef.ok<-!is.na(gamma.c)
                psi<-psi[id.coef.ok]
                if(length(psi)<=0) {
                    warning(paste("All breakpoints have been removed after",it,"iterations.. returning 0"), call. = FALSE)
                    return(0)
                }
                gamma.c<-gamma.c[id.coef.ok]
                beta.c<-beta.c[id.coef.ok]
                Z<-Z[, id.coef.ok, drop=FALSE]
                rangeZ <- rangeZ[,id.coef.ok, drop=FALSE]
                limZ <- limZ[,id.coef.ok, drop=FALSE]
                nomiOK<-nomiOK[id.coef.ok] #salva i nomi delle U per i psi ammissibili
                id.psi.group<-id.psi.group[id.coef.ok]
                names(psi)<-id.psi.group
            }
        }
        psi.old<-psi
        psi <- psi.old + gamma.c/beta.c
        if(!is.null(digits)) psi<-round(psi, digits)
        PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
        
        #--modello con il nuovo psi
        U<-(Z-PSI)*(Z>PSI) #in seg.(g)lm.fit l'ho chiamata U..
        if(pow[1]!=1) U<-U^pow[1]
        #call.noV$init<-quote(coef(...........))
        obj1 <- suppressWarnings(try(eval(call.noV), silent=TRUE))
        if(class(obj1)[1]=="try-error")  obj1 <- suppressWarnings(eval(call.noV.noinit))
        L1<- if(class(obj1)[1]=="try-error") L0+10 else -obj1$loglik
        use.k<-k<-1
        L1.k<-NULL
        L1.k[length(L1.k)+1]<-L1
        
        while(L1>L0){
            k<-k+1
            use.k <- if(useExp.k) 2^(k-1) else k
            #     if(k>=4){
            #        xx<-1:k
            #        use.k<-est.k(xx, -L1.k[1:k],-L0)
            #      }
            psi <- psi.old + (gamma.c/beta.c)/(use.k*h)
            #psi <- psi.old[id.psi.ok] + (gamma.c[id.psi.ok]/beta.c[id.psi.ok])/(use.k*h)
            if(!is.null(digits)) psi<-round(psi, digits)
            PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
            #qui o si aggiusta psi per farlo rientrare nei limiti, o si elimina, oppure se obj1 sotto non funziona semplicemente continua..
            U<-(Z-PSI)*(Z>PSI)
            if(pow[1]!=1) U<-U^pow[1]
            obj1 <- suppressWarnings(try(eval(call.noV), silent=TRUE))
            L1<- if(class(obj1)[1]=="try-error") L0+10 else -obj1$loglik
            L1.k[length(L1.k)+1]<-L1
            if(1/(use.k*h)<min.step){
                #        #        #warning("step halving too small") 
                break}
        } #end while L0-L1
        
        if (visual) {
            flush.console()
            cat(paste("iter = ", sprintf("%2.0f",it),
                      "  llik = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), -L1), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
                      "  k = ", sprintf("%2.0f", k),
                      "  n.psi = ",formatC(length(unlist(psi)),digits=0,format="f"), 
                      "  est.psi = ",paste(formatC(unlist(psi),digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                      sep=""), "\n")
        }
        
        epsilon <- if(conv.psi) max(abs((psi -psi.old)/psi.old)) else (L0 - L1)/(abs(L0) + 0.1) 
        L0<-L1
        #U <-U1 #le ho sempre chiamate U, non c'e' bisogno di ri-assegnare..
        
        k.values[length(k.values)+1]<-use.k
        psi.values[[length(psi.values) + 1]] <- psi
        dev.values[length(dev.values) + 1] <- L0
        
        #Mi sa che non servono i controlli.. soprattutto se non ha fatto step-halving
        #check if i psi ottenuti sono nel range o abbastanza lontani
        id.psi.far <-far.psi(Z, PSI, id.psi.group, TRUE, fc=opz$fc)
        id.psi.in <- in.psi(limZ, PSI, TRUE)
        id.psi.ok <- id.psi.in & id.psi.far 
        
        if(!all(id.psi.ok)){
            if(fix.npsi){
#                psi<-adj.psi(psi, limZ) #within range!!!
#                id.psi.far<-far.psi(Z, PSI, id.psi.group, TRUE)
                psi<-psi*ifelse(id.psi.far,1,.9) 
                PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
                id.psi.changed[it]<-TRUE
                    } else {
                Z<-Z[, id.psi.ok, drop=FALSE]
                PSI<-PSI[, id.psi.ok, drop=FALSE]
                rangeZ <- rangeZ[,id.psi.ok,drop=FALSE]
                limZ <- limZ[,id.psi.ok,drop=FALSE]
                nomiOK<-nomiOK[id.psi.ok] #salva i nomi delle U per i psi ammissibili
                id.psi.group<-id.psi.group[id.psi.ok]
                psi.old<- psi.old[id.psi.ok]
                psi<- psi[id.psi.ok]  
                names(psi)<-id.psi.group
                if(ncol(PSI)<=0) {
                    warning(paste("All breakpoints have been removed after",it,"iterations.. returning 0"), call. = FALSE)
                    return(0)
                }
            }
        }
        if (it >= it.max) {
            id.warn <- TRUE
            break
        }
    } #end while_it
    ##=============================================================================
    if(id.psi.changed[length(id.psi.changed)]) warning(paste("Some psi (", (1:length(psi))[!id.psi.far],
                                                             ") changed after the last iter.",sep=""), call. = FALSE)
    if(id.warn) warning(paste("max number of iterations (", it,") attained",sep=""), call. = FALSE)
    
    attr( psi.values, "dev") <- dev.values
    attr( psi.values, "k")<- k.values
    
    #ordina i breakpoints.. 
    psi<-unlist(tapply(psi, id.psi.group, sort))
    names(psi)<-id.psi.group

    names.coef<-names(coef(obj)) #names(obj$coefficients) #obj e' quello vecchio che include U1,.. V1,...
    
    PSI.old<-PSI
    PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
    
    #U e V possono essere cambiati (rimozione/ordinamento psi.. ) per cui si deve ricalcolare il tutto, altrimenti sarebbe uguale a U1 e obj1
    if(sd(PSI-PSI.old)>0 || id.psi.changed[length(id.psi.changed)]){
        U <- (Z-PSI)*(Z>PSI)
        colnames(U)<-paste("U", 1:ncol(U), sep = "")
        V <- -(Z>PSI)
        colnames(V)<-paste("V", 1:ncol(V), sep = "")
        #X <- cbind(XREG, U, V)
        #rownames(X) <- NULL
        obj <- suppressWarnings(try(eval(call.noV), silent=TRUE))
        L1<- -obj$loglik
    } else {
        obj<-obj1
    }
    obj$coef<-c(obj$coef, rep(0,ncol(V)))
    names(obj$coef)<-names.coef
    obj$epsilon <- epsilon
    obj$it <- it
    
    obj<-list(obj=obj,it=it,psi=psi, psi.values=psi.values, U=U,V=V,rangeZ=rangeZ,
              epsilon=epsilon,nomiOK=nomiOK, SumSquares.no.gap=L1, id.psi.group=id.psi.group,id.warn=id.warn) #inserire id.psi.ok?
    return(obj)
}

# if(return.all.sol) {
        #     obj.noV <- suppressWarnings(eval(call.noV)) #, envir=mfExt
        #     #mio.init.noV<-obj.noV$coef
        #     #mio.init.noV<- c(0,obj.noV$coef[-1])
        #     dev.new1 <- -obj.noV$loglik 
        #     #dev.new1 <- sum(mylm(x = cbind(XREG, U), y = y, w = w, offs = offs)$residuals^2)
        #     }
        # 
        # beta.c<-coef(obj)[nomiU]
        # gamma.c<-coef(obj)[nomiV]

