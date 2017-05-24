seg.def.fit<-function(obj, Z, PSI, mfExt, opz, return.all.sol=FALSE){
#-----------------
fn.costr<-function(n.psi,isLeft=1,isInterc=1){
#build the constraint matrix
#isLeft: TRUE (or 1) if there is the left slope
#isInterc: TRUE (or 1) if there is the intercept..
IU<- -diag(n.psi)
sumU<- diag(n.psi) #n. of diff slopes
sumU[row(sumU)>col(sumU)]<-1
if(isLeft) {
    sumU<-cbind(1, sumU)
    IU<-diag(c(1, -rep(1, n.psi)))
    }
    A<-rbind(IU,sumU)
    if(isInterc) {
      A<-rbind(0,A)
      A<-cbind(c(1, rep(0,nrow(A)-1)), A)
      } 
    #add zeros for the V coeffs
    A<-cbind(A, matrix(0,nrow(A), n.psi))
    A
    }
#-----------------
dpmax<-function(x,y,pow=1){
#deriv pmax
        if(pow==1) ifelse(x>y, -1, 0)
         else -pow*pmax(x-y,0)^(pow-1)
         }
#-----------
    vincoli<- FALSE
    c1 <- apply((Z <= PSI), 2, all)
    c2 <- apply((Z >= PSI), 2, all)
    if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2))) stop("psi out of the range")
    #
    digits<-opz$digits
    pow<-opz$pow
    nomiOK<-opz$nomiOK
    toll<-opz$toll
    h<-opz$h
    gap<-opz$gap
    stop.if.error<-opz$stop.if.error
    dev.new<-opz$dev0
    visual<-opz$visual
    id.psi.group<-opz$id.psi.group
    it.max<-old.it.max<-opz$it.max
    rangeZ <- apply(Z, 2, range)
    psi<-PSI[1,]
    names(psi)<-id.psi.group
    #H<-1
    it <- 1
    epsilon <- 10
    dev.values<-psi.values <- NULL
    id.psi.ok<-rep(TRUE, length(psi))
    
    nomiU<- opz$nomiU
    nomiV<- opz$nomiV
    call.ok <- opz$call.ok
    call.noV <- opz$call.noV

#browser()
    if(is.null(opz$constr)) opz$constr<-0
    if((opz$constr %in% 1:2) && class(obj)=="rq"){
      vincoli<-TRUE
      call.ok$method<-"fnc"
      call.ok$R<-quote(R)
      call.ok$r<-quote(r)

      call.noV$method<-"fnc"
      call.noV$R<-quote(R.noV)
      call.noV$r<-quote(r)
      }
    
    fn.obj<-opz$fn.obj
    toll<-opz$toll
    k<-ncol(Z)
    while (abs(epsilon) > toll) {
        #k<-ncol(Z)
        U <- pmax((Z - PSI), 0)^pow[1]#U <- pmax((Z - PSI), 0)
        V <- dpmax(Z,PSI,pow=pow[2])# ifelse((Z > PSI), -1, 0)
        for(i in 1:k) {
          mfExt[nomiU[i]] <- U[,i]
          mfExt[nomiV[i]] <- V[,i]
        }
        R<- fn.costr(ncol(U),1,1)
        R.noV<-R[,-((ncol(R)-1)+seq_len(ncol(U))),drop=FALSE]
        r<- rep(0, nrow(R))
        obj <- suppressWarnings(eval(call.ok, envir=mfExt))
        dev.old<-dev.new
        dev.new <- dev.new1 <- eval(parse(text=fn.obj), list(x=obj)) #control$f.obj should be something like "sum(x$residuals^2)" or "x$dev"   
        if(length(dev.new)<=0) stop("error in the objective to be minimized, see 'fn.obj'")
        if(return.all.sol) {
            obj.noV <- suppressWarnings(eval(call.noV, envir=mfExt))
            dev.new1 <- eval(parse(text=fn.obj), list(x=obj.noV))
            #dev.new1 <- sum(mylm(x = cbind(XREG, U), y = y, w = w, offs = offs)$residuals^2)
            }
        dev.values[[length(dev.values) + 1]] <- dev.new1
        if (visual) {
            flush.console()
            if (it == 1)
                cat(0, " ", formatC(dev.old, 3, format = "f"),
                  "", "(No breakpoint(s))", "\n")
            spp <- if (it < 10) "" else NULL
            cat(it, spp, "", formatC(dev.new, 3, format = "f"), "",length(psi),"\n")
            #cat(paste("iter = ", it, spp," dev = ",formatC(dev.new,digits=3,format="f"), " n.psi = ",formatC(length(psi),digits=0,format="f"), sep=""), "\n")
        }
        epsilon <- (dev.new - dev.old)/(dev.old + .001)
        obj$epsilon <- epsilon
        it <- it + 1
        obj$it <- it
        beta.c<-coef(obj)[nomiU]
        gamma.c<-coef(obj)[nomiV]
        
        if (it > it.max) break
        psi.values[[length(psi.values) + 1]] <- psi.old <- psi
 #       if(it>=old.it.max && h<1) H<-h
        psi <- psi.old + h*gamma.c/beta.c
        if(!is.null(digits)) 
        psi<-round(psi, digits)
        PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
        #check if psi is admissible..
        a <- apply((Z <= PSI), 2, all) #prima era solo <
        b <- apply((Z >= PSI), 2, all) #prima era solo >
        if(stop.if.error) {
          isErr<- (sum(a + b) != 0 || is.na(sum(a + b)))
            if(isErr) {
              if(return.all.sol) return(list(dev.values, psi.values)) else stop("(Some) estimated psi out of its range")
              }
            } else {
          id.psi.ok<-!is.na((a+b)<=0)&(a+b)<=0
          Z <- Z[,id.psi.ok,drop=FALSE]
          psi <- psi[id.psi.ok]
          PSI <- PSI[,id.psi.ok,drop=FALSE]
          
          ToDeletenomiU<-nomiU[!id.psi.ok] #salva i nomi delle U per i psi ammissibili
          ToDeletenomiV<-nomiV[!id.psi.ok] #salva i nomi delle V per i psi ammissibili
          if(length(ToDeletenomiU)>0 || length(ToDeletenomiV)>0) {for(nn in c(ToDeletenomiU, ToDeletenomiV)) {mfExt[[nn]]<-NULL}}

          nomiOK<-nomiOK[id.psi.ok] #salva i nomi delle U per i psi ammissibili
          nomiU<-nomiU[id.psi.ok] #salva i nomi delle U per i psi ammissibili
          nomiV<-nomiV[id.psi.ok] #salva i nomi delle V per i psi ammissibili
                    
          id.psi.group<-id.psi.group[id.psi.ok]
          names(psi)<-id.psi.group
          if(ncol(PSI)<=0) return(0)
          k<-ncol(Z)
          #aggiorna la call, altrimenti il modello avra' sempre lo stesso numero di termini anche se alcuni psi vengono rimossi!!!
          Fo <- update.formula(opz$formula.orig, as.formula(paste(".~.+", paste(c(nomiU, nomiV), collapse = "+"))))
          Fo.noV <- update.formula(opz$formula.orig, as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
   
          
          call.ok <- update(obj, formula = Fo,  evaluate=FALSE, data = mfExt) 
          call.noV <- update(obj, formula = Fo.noV,  evaluate=FALSE, data = mfExt) 

        
        } #end else
        #obj$psi <- psi
    } #end while

    psi<-unlist(tapply(psi, id.psi.group, sort))
    names(psi)<-id.psi.group
    PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
    #aggiunto da qua..
    

    U <- pmax((Z - PSI), 0)
    V <- ifelse((Z > PSI), -1, 0)
    for(i in 1:k) {
          mfExt[nomiU[i]] <- U[,i]
          mfExt[nomiV[i]] <- V[,i]
          }

##LA DOMANDA E': PERCHE' QUI STIMA UN MODELLO SENZA V SE POI VIENE RISTIMATO in segmented.default (o segmented.lm o segmented.glm?)
##RE: il valore di SS.new serve per il boot restart.
#Invece la domanda e': non si puo' restituire direttamente obj.new senza bisogno di sostituire i valori in obj ?
    obj.new <- suppressWarnings(eval(call.noV, envir=mfExt))
    SS.new <- eval(parse(text=fn.obj), list(x=obj.new)) #sum(obj.new$residuals^2)
    if(!gap){
          obj<-obj.new
          #names.coef<-names(obj$coefficients)
          #obj$coefficients<-c(obj.new$coefficients, rep(0,ncol(V)))
          #names(obj$coefficients)<-names.coef
          #obj$residuals<-obj.new$residuals
          #obj$fitted.values<-obj.new$fitted.values
          #obj$linear.predictors<-obj.new$linear.predictors
          #obj$deviance<-obj.new$deviance
          #obj$weights<-obj.new$weights
          #obj$aic<-obj.new$aic #+ 2*ncol(V) #ho fatto la modifica in segmented.glm(): "objF$aic<-obj$aic + 2*k"
          } else {
          obj <- suppressWarnings(eval(call.ok, envir=mfExt))
          }
    obj$epsilon <- epsilon
    obj$it <- it
    #fino a qua..
    obj<-list(obj=obj,it=it,psi=psi,psi.values=psi.values,U=U,V=V,rangeZ=rangeZ,
        epsilon=epsilon,nomiOK=nomiOK, SumSquares.no.gap=SS.new, id.psi.group=id.psi.group, 
        nomiV=nomiV, nomiU=nomiU, mfExt=mfExt) #inserire id.psi.ok?
    #browser()
    if(vincoli) {
        obj$R<-R
        obj$R.noV<-R.noV
        obj$r<-r
        } 
     return(obj)
    }

