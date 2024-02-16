seg.Ar.fit.boot<-function(obj, XREG, Z, PSI, opz, n.boot=10, size.boot=NULL, jt=FALSE,
    nonParam=TRUE, random=FALSE, break.boot=n.boot){
#random se TRUE prende valori random quando e' errore: comunque devi modificare qualcosa (magari con it.max)
#     per fare restituire la dev in corrispondenza del punto psi-random
#nonParm. se TRUE implemneta il case resampling. Quello semiparam dipende dal non-errore di
      extract.psi<-function(lista){
      #serve per estrarre il miglior psi..
        dev.values<-lista[[1]][-1] #remove the 1st one referring to model without psi
        psi.values<-lista[[2]][-1] #remove the 1st one (NA)
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
        mY <- mean(obj$residuals)
        vv <- strsplit(paste(strsplit(paste(mY),"\\.")[[1]], collapse=""),"")[[1]]
        vv<-vv[vv!="0"]
        vv=na.omit(vv[1:5])
        seed <-eval(parse(text=paste(vv, collapse="")))
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
      opz.boot<-opz
      opz1<-opz
      
      opz.boot$pow=c(1,1) #c(1.1,1.2)
      opz.boot$it.max<-20
      
      opz1$it.max <-0
      n<-nrow(Z)
      o0<-try(suppressWarnings(seg.Ar.fit(obj, XREG, Z, PSI, opz)), silent=TRUE)
      rangeZ <- apply(Z, 2, range) #serve sempre
      alpha <- opz$alpha
      limZ <- apply(Z, 2, quantile, names = FALSE, probs = c(alpha, 1 - alpha))
      
      if(!is.list(o0)) {
          o0<- seg.Ar.fit(obj, XREG, Z, PSI, opz, return.all.sol=TRUE)
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
            PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
            o0<-try(suppressWarnings(seg.Ar.fit(obj, Z, PSI1, opz1)), silent=TRUE)
            ss00<-o0$SumSquares.no.gap
          } else {
          est.psi00<-est.psi0<-apply(PSI,2,mean)
          ss00<-opz$dev0
        }
        }

      n.intDev0<-nchar(strsplit(as.character(ss00),"\\.")[[1]][1])
      
      all.est.psi.boot<-all.selected.psi<-all.est.psi<-matrix(, nrow=n.boot, ncol=length(est.psi0))
      all.ss<-all.selected.ss<-rep(NA, n.boot)
      if(is.null(size.boot)) size.boot<-n

      Z.orig<-Z
      #if(visualBoot) cat(0, " ", formatC(opz$dev0, 3, format = "f"),"", "(No breakpoint(s))", "\n")
      count.random<-0
      for(k in seq(n.boot)){
        ##se gli *ultimi* n.boot.rev valori di ss sono uguali, cambia i psi...
        n.boot.rev<- 3 #3 o 4?
        diff.selected.ss <- rev(diff(na.omit(all.selected.ss)))
        #if(length(na.omit(diff(all.selected.ss[1:n.boot.rev])))==(n.boot.rev-1) && all(round(diff(all.selected.ss[1:n.boot.rev]),6)==0)){
        if(length(diff.selected.ss)>=(n.boot.rev-1) && all(round(diff.selected.ss[1:(n.boot.rev-1)],6)==0)){
          qpsi<-sapply(1:ncol(Z),function(i)mean(est.psi0[i]>=Z[,i]))
          qpsi<-ifelse(abs(qpsi-.5)<.1,.8,qpsi)
          est.psi0<-sapply(1:ncol(Z),function(i)quantile(Z[,i],probs=1-qpsi[i],names=FALSE))
        }
        
        PSI <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
        if(jt) Z<-apply(Z.orig,2,jitter)
        if(nonParam){
            id<-sample(n, size=size.boot, replace=TRUE)
            o.boot<-try(suppressWarnings(seg.Ar.fit(obj, XREG[id,,drop=FALSE], Z[id,,drop=FALSE], PSI[id,,drop=FALSE], opz.boot)), silent=TRUE)
          } else {
            yy<-fitted.ok+sample(residuals(o0), size=n, replace=TRUE)
##---->     o.boot<-try(seg.lm.fit(yy, XREG, Z.orig, PSI, weights, offs, opz.boot), silent=TRUE)
            #in realta' la risposta dovrebbe essere "yy" da cambiare in mfExt
            o.boot<- try(suppressWarnings(seg.Ar.fit(obj, XREG, Z.orig, PSI, opz.boot)), silent=TRUE)
          }
          if(is.list(o.boot)){
            all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi
            } else {
            est.psi.boot<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
            }
          PSI <- matrix(rep(est.psi.boot, rep(nrow(Z), length(est.psi.boot))), ncol = length(est.psi.boot))
          #opz$h<-max(opz$h*.9, .2)
          opz$it.max<-opz$it.max+1
          o <- try(suppressWarnings(seg.Ar.fit(obj, XREG, Z.orig, PSI, opz, return.all.sol=TRUE)), silent=TRUE)
          if(!is.list(o) && random){
              est.psi0<-apply(limZ,2,function(r) runif(1,r[1],r[2]))
              PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
              o <- try(suppressWarnings(seg.Ar.fit(obj, XREG, Z, PSI1, opz1)), silent=TRUE)
              count.random<-count.random+1
              }
          if(is.list(o)){
              if(!"coef"%in%names(o$obj)) o<-extract.psi(o)
              all.est.psi[k,]<-o$psi
              all.ss[k]<-o$SumSquares.no.gap
              if(o$SumSquares.no.gap<=ifelse(is.list(o0), o0$SumSquares.no.gap, 10^12)) o0<-o
              est.psi0<-o0$psi
              all.selected.psi[k,] <- est.psi0
              all.selected.ss[k]<-o0$SumSquares.no.gap #min(c(o$SumSquares.no.gap, o0$SumSquares.no.gap))
              }
          if (visualBoot) {
              flush.console()
              #n.intDev0<-nchar(strsplit(as.character(dev.values[2]),"\\.")[[1]][1])
              cat(paste("boot sample = ", sprintf("%2.0f",k),
                        "  opt.llik = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), -o0$SumSquares.no.gap), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
                        "  n.psi = ",formatC(length(unlist(est.psi0)),digits=0,format="f"), 
                        "  est.psi = ",paste(formatC(unlist(est.psi0),digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                        sep=""), "\n")
                  }
          asss<-na.omit(all.selected.ss)
          if(length(asss)>break.boot){
            if(all(rev(round(diff(asss),6))[1:(break.boot-1)]==0)) break
          }
        } #end n.boot
      all.selected.psi<-rbind(est.psi00,all.selected.psi)
      all.selected.ss<-c(ss00, all.selected.ss)

      ris<-list(all.selected.psi=drop(all.selected.psi),all.selected.ss=all.selected.ss, all.psi=all.est.psi, all.ss=all.ss)

      if(is.null(o0$obj)){
          PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
          o0 <- try(suppressWarnings(seg.Ar.fit(obj, XREG, Z, PSI1, opz1)), silent=TRUE)
          warning("The final fit can be unreliable (possibly mispecified segmented relationship)", call.=FALSE, immediate.=TRUE)
      }
      if(!is.list(o0)) return(0)
      o0$boot.restart<-ris
      o0$seed<-seed
      #rm(.Random.seed, envir=globalenv())
      return(o0)
      }