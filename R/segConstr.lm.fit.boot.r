segConstr.lm.fit.boot <- function(y, XREG, Z, PSI, w, offs, opz, n.boot=10, size.boot=NULL, jt=FALSE,
    nonParam=TRUE, random=FALSE, break.boot=n.boot){
#random se TRUE prende valori random quando e' errore: comunque devi modificare qualcosa (magari con it.max)
#     per fare restituire la dev in corrispondenza del punto psi-random
#nonParm. se TRUE implemneta il case resampling. Quello semiparam dipende dal non-errore di
#----------------------------------
#  sum.of.squares<-function(obj.seg){
#      #computes the "correct" SumOfSquares from a segmented" fit
#      b<-obj.seg$obj$coef
#      X<-qr.X(obj.seg$obj$qr) #X<-model.matrix(obj.seg)
#      X<-X[,!is.na(b)]
#      b<-b[!is.na(b)]
#      rev.b<-rev(b)
#      rev.b[1:length(obj.seg$psi)]<-0
#      b<-rev(rev.b)
#      new.fitted<-drop(X%*%b)
#      new.res<- obj.seg$obj$residuals + obj.seg$obj$fitted - new.fitted
#      ss<-sum(new.res^2)
#      ss
#      }
#--------
  extract.psi<-function(lista){
    #serve per estrarre il miglior psi..
    dev.values<-lista[[1]]
    psi.values<-lista[[2]]
    if(any(is.na(psi.values[[1]]))) {#se la 1 componente e' NA (fino alla versione 2.0-3 era cosi'... perche' in dev.values c'erano 
      #  anche i valori relativi al modello senza psi.. )
      dev.values<-dev.values[-1] #remove the 1st one referring to model without psi
      psi.values<-psi.values[-1]
    }
    dev.ok<-min(dev.values)
    id.dev.ok<-which.min(dev.values)
    if(is.list(psi.values))  psi.values<-matrix(unlist(psi.values),
                                                nrow=length(dev.values), byrow=TRUE)
    if(!is.matrix(psi.values)) psi.values<-matrix(psi.values)
    psi.ok<-psi.values[id.dev.ok,]
    r<-list(SumSquares.no.gap=dev.ok, psi=psi.ok)
    r
  }
  #-------------
      if(is.null(opz$seed)){
        mY <- mean(y)
        sepDec<-if(options()$OutDec==".") "\\." else "\\,"
        vv <- strsplit(paste(strsplit(paste(mY), sepDec)[[1]], collapse=""),"")[[1]]
        vv<-vv[vv!="0"]
        vv=na.omit(vv[1:5])
        seed <-eval(parse(text=paste(vv, collapse="")))
        if(is.null(seed)) seed <- 1
        set.seed(seed)
      } else {
        if(is.na(opz$seed)) {
          seed <-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
          set.seed(seed)
        } else {
          seed <-opz$seed
          set.seed(opz$seed)
        }
      }  
      
      visualBoot<-opz$visualBoot
      #opz.boot<-opz
      #opz.boot$pow=c(1,1) #c(1.1,1.2)
      opz1<-opz
      opz1$it.max <- 0
      n<-length(y)
      rangeZ <- apply(Z, 2, range) #serve sempre
      
      alpha <- opz$alpha
      #limZ <- apply(Z, 2, quantile, names = FALSE, probs = c(alpha[1], alpha[2]))
      limZ <- if(is.null(opz$limZ)) apply(Z, 2, quantile, names=FALSE, probs=c(alpha[1],alpha[2])) else opz$limZ
      
      o0<-try(suppressWarnings(segConstr.lm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=FALSE)), silent=TRUE)
      if(!is.list(o0)){
        o0<-try(suppressWarnings(segConstr.lm.fit(y, XREG, Z, opz$PSI1, w, offs, opz, return.all.sol=FALSE)), silent=TRUE)
      }
      
      if(!is.list(o0)) {
          o0<- suppressWarnings(segConstr.lm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=TRUE))
          o0<-extract.psi(o0)
          ss00<-opz$dev0
          if(!nonParam) {warning("using nonparametric boot");nonParam<-TRUE}
          }
      if(is.list(o0)){
        est.psi00<-est.psi0<-o0$psi
        ss00<-o0$SumSquares.no.gap
        if(!nonParam) fitted.ok<-fitted(o0)
        } else {
          if(!nonParam) stop("the first fit failed and I cannot extract fitted values for the semipar boot")
          if(random) {
            est.psi00<-est.psi0<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
            PSI1 <- matrix(est.psi0, n, ncol = length(est.psi0), byrow = TRUE)
            o0<-try(suppressWarnings(segConstr.lm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
            ss00<-o0$SumSquares.no.gap
          } else {
          est.psi00<-est.psi0<-apply(PSI,2,mean)
          ss00<-opz$dev0
        }
        }

      n.intDev0<-nchar(strsplit(as.character(ss00),options()$OutDec)[[1]][1])
      
      all.est.psi.boot<-all.selected.psi<-all.est.psi<-matrix(NA, nrow=n.boot, ncol=length(est.psi0))
      all.ss<-all.selected.ss<-rep(NA, n.boot)
      if(is.null(size.boot)) size.boot<-n
#      na<- ,,apply(...,2,function(x)mean(is.na(x)))
      Z.orig<-Z
#     if(visualBoot) cat(0, " ", formatC(opz$dev0, 3, format = "f"),"", "(No breakpoint(s))", "\n")
      count.random<-0
      id.uguali<-0
      k.psi.change<- 1
      alpha<-.1
      n.boot.rev<- 3 #3 o 4?
      
      for(k in seq(n.boot)){
        #if(k==4) browser()
        ##se gli *ultimi* n.boot.rev valori di ss sono uguali, cambia i psi...
        diff.selected.ss <- rev(diff(na.omit(all.selected.ss)))
        #if(length(na.omit(diff(all.selected.ss[1:n.boot.rev])))==(n.boot.rev-1) && all(round(diff(all.selected.ss[1:n.boot.rev]),6)==0)){
        if(length(diff.selected.ss)>=(n.boot.rev-1) && all(round(diff.selected.ss[1:(n.boot.rev-1)],6)==0)){
          qpsi<-sapply(1:ncol(Z),function(i)mean(est.psi0[i]>=Z[,i]))
          qpsi<-ifelse(abs(qpsi-.5)<.1, alpha, qpsi)
          alpha<-1-alpha
          est.psi0<-sapply(1:ncol(Z),function(i)quantile(Z[,i],probs=1-qpsi[i],names=FALSE))
        }
        ########################### 25/7/24 #####
        est.psi0 <- unlist(tapply(est.psi0, opz$id.psi.group, sort))
        #########################################
        PSI <- matrix(est.psi0, n, ncol = length(est.psi0), byrow=TRUE)
        if(jt) Z<-apply(Z.orig,2,jitter)
        if(nonParam){
              id<-sample(n, size=size.boot, replace=TRUE)
              o.boot<-try(suppressWarnings(segConstr.lm.fit(y[id], XREG[id,,drop=FALSE], Z[id,,drop=FALSE], PSI[id,,drop=FALSE],
                w[id], offs[id], opz)), silent=TRUE)
        } else {
              yy<-fitted.ok+sample(residuals(o0),size=n, replace=TRUE)
              o.boot<-try(suppressWarnings(segConstr.lm.fit(yy, XREG, Z.orig, PSI, weights, offs, opz)), silent=TRUE)
        }
        if(is.list(o.boot)){
            all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi
        } else {
            est.psi.boot<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
            est.psi.boot<- unlist(tapply(est.psi.boot, opz$id.psi.group, sort))
        }
        #if(k==7) browser()
        ### se est.psi.boot non e' cambiato (e puoi vederlo da all.est.psi.boot), allora cambialo!
        
        
        PSI <- matrix(est.psi.boot, n, ncol = length(est.psi.boot), byrow=TRUE)
        #opz$h<-max(opz$h*.9, .2)
        opz$it.max<-opz$it.max+1
        o<-try(suppressWarnings(segConstr.lm.fit(y, XREG, Z.orig, PSI, w, offs, opz, return.all.sol=TRUE)), silent=TRUE)
        if(!is.list(o) && random){
                est.psi0<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
                PSI1 <- matrix(est.psi0, n, ncol = length(est.psi0), byrow=TRUE)
                o<-try(suppressWarnings(segConstr.lm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
                count.random<-count.random+1
        }
        #se il modello e' stato stimato controlla se la soluzione e' migliore..
        if(is.list(o)){
              if(!"coefficients"%in%names(o$obj)) o<-extract.psi(o)
              all.est.psi[k,]<-o$psi
              all.ss[k]<-o$SumSquares.no.gap
              if(o$SumSquares.no.gap<=ifelse(is.list(o0), o0$SumSquares.no.gap, 10^12)) {o0<-o; k.psi.change<- k}
              est.psi0<-o0$psi
              all.selected.psi[k,] <- est.psi0
              all.selected.ss[k]<-o0$SumSquares.no.gap #min(c(o$SumSquares.no.gap, o0$SumSquares.no.gap))
        }
            
        if(visualBoot) {
              flush.console()
              #      spp <- if (it < 10) " " else NULL
              #      cat(paste("iter = ", spp, it,
              #                "  dev = ",sprintf('%8.5f',L1), #formatC(L1,width=8, digits=5,format="f"), #era format="fg"
              #n.intDev0<-nchar(strsplit(as.character(dev.values[2]),"\\.")[[1]][1])
              cat(paste("boot sample = ", sprintf("%2.0f",k),
                    #"  opt.dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), o0$SumSquares.no.gap), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
                    "  opt.dev = ", sprintf("%1.5f", as.numeric(strsplit(format(o0$SumSquares.no.gap, scientific=TRUE), "e")[[1]][1])),
                    "  n.psi = ",formatC(length(unlist(est.psi0)),digits=0,format="f"), 
                    "  est.psi = ",paste(formatC(unlist(est.psi0),digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                    sep=""), "\n")
        }
        #conta i valori ss uguali.. cosi puoi fermarti prima..
        asss<-na.omit(all.selected.ss)
        if(length(asss)>break.boot){
          if(all(rev(round(diff(asss),6))[1:(break.boot-1)]==0)) break
        }
        #id.uguali<-(round(diff(all.selected.ss[c(k-1,k-2)]),6)==0)+id.uguali      
        } #end n.boot

      

      all.selected.psi<-rbind(est.psi00,all.selected.psi)
      all.selected.ss<-c(ss00, all.selected.ss)


      # SS.ok<-min(all.selected.ss)
      # id.accept<- ((abs(all.ss-SS.ok)/SS.ok )<= 0.05)
      # psi.mean<-apply(all.est.psi[id.accept,,drop=FALSE], 2, mean)
      # est.psi0<-psi.mean
      # devi ristimare il modello con psi.mean
      # PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
      # o0<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)

      ris<-list(all.selected.psi=drop(all.selected.psi),all.selected.ss=all.selected.ss, all.psi=all.est.psi, all.ss=all.ss)

      if(is.null(o0$obj)){
        #quando vengono restituiti psi troppo vicini e l'SE non si puo' calcolare, possiamo distanziarli..
        #Pero' il processo deve essere esteso nel caso in cui ci sono 3 psi vicini..
        min.n <- opz$min.n-1
        if(min.n>1){
          min1<- function(x, k=min.n-1){
            for(i in 1:k) x<-x[-which.min(x)]
            min(x)
          }
          max1<-function(x,k=min.n-1){
            for(i in 1:k) x<-x[-which.max(x)]
            max(x)
          }
        } else {
          min1<-min
          max1<-max
        }
        npsi <- tapply(opz$id.psi.group, opz$id.psi.group, length)
        nomiAll <- colnames(rangeZ) #rep(opz$nomiSeg, npsi)
        nomiSeg <- unique(nomiAll)
        newPsi<-vector("list", length(npsi) )
        for(.j in 1:length(npsi)){
          psi.j <- sort(est.psi0[opz$id.psi.group==.j]) #psi della stessa variabile segmented
          id  <- nomiSeg[.j]==nomiAll
          Z.ok <- unique(Z[, id, drop=FALSE][,1])
          m.j <- min(limZ[1,id])
          M.j <- max(limZ[2,id])
          #h=1/1.05
          for(.k in 1:length(psi.j)){
            id.group<-cut(Z.ok, c(m.j-10^8, psi.j, M.j+10^8), labels=FALSE)
            n.j<-tabulate(id.group)#<=min.n
            #per ogni psi calcola il min e il max dei segmenti prima e dopo psi. 
            #se questi segmenti hanno min.n osservazioni considera u min e max fittizzi per evitare che il nuovo psi
            #modificato porti a segmenti con bassa numerosita'..
            M.j.k<- if(n.j[.k]>0) max1(Z.ok[id.group==.k])  -10^6*(n.j[.k]<=min.n) else -10^6*(n.j[.k]<=min.n) 
            m.j.k<- if(n.j[.k+1]>0) min1(Z.ok[id.group==.k+1])+10^6*(n.j[.k+1]<=min.n) else  10^6*(n.j[.k]<=min.n)
            psi.j[.k]<- psi.j[.k] + ifelse(abs(M.j.k-psi.j[.k])<abs(m.j.k-psi.j[.k]), M.j.k-psi.j[.k]-.0001, m.j.k-psi.j[.k]+.0001  )
          }
          newPsi[[.j]]<-psi.j
        } #end .j
        est.psi0 <- unlist(newPsi)
        PSI1 <- matrix(est.psi0, n, ncol = length(est.psi0), byrow=TRUE)
        o0<-try(suppressWarnings(segConstr.lm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
        warning("Breakpoint estimates have been outdistanced to allow finite estimates and st.errs", call.=FALSE, immediate.=TRUE)
        #warning(" 'The final fit (if returned) could be unreliable. Reduce no. of psi or try to increase 'break.boot'", call.=FALSE, immediate.=TRUE)
        #warning("'Convergence' is suspect: the final fit could be unreliable. Try to re-run by increasing 'break.boot'", call.=FALSE, immediate.=TRUE)
      }
      
      
      
      if(!is.list(o0)) return(0)
      o0$boot.restart<-ris
      o0$seed <- seed
      #rm(.Random.seed, envir=globalenv())
      return(o0)
      }