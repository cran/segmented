segConstr.glm.fit.boot <- function(y, XREG, Z, PSI, w, offs, opz, n.boot=10, size.boot=NULL, jt=FALSE,
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
    	dev.values<-lista[[1]][-1] #remove the 1st one referring to model without psi
    	psi.values<-lista[[2]][-1] #remove the 1st one (NA)
    	dev.ok<-min(dev.values)
    	id.dev.ok<-which.min(dev.values)
    	if(is.list(psi.values))  psi.values<-matrix(unlist(psi.values),
    		nrow=length(dev.values), byrow=TRUE)
    	if(!is.matrix(psi.values)) psi.values<-matrix(psi.values)
    	psi.ok<-psi.values[id.dev.ok,]
    	r<-list(dev.no.gap=dev.ok, psi=psi.ok)
    	r
	}
#-------------
      visualBoot<-opz$visualBoot
      opz.boot<-opz
      opz.boot$pow=c(1,1) #c(1.1,1.2)
      opz1<-opz
      opz1$it.max <-1
      n<-length(y)
      o0<-try(suppressWarnings(segConstr.glm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=FALSE)), silent=TRUE)
      rangeZ <- apply(Z, 2, range) #serve sempre
      
      alpha <- opz$alpha
      limZ <- apply(Z, 2, quantile, names = FALSE, probs = c(alpha[1], alpha[2]))
      
      if(!is.list(o0)) {
          o0<- suppressWarnings(segConstr.glm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=TRUE))
          o0<-extract.psi(o0)
          ss00<-opz$dev0
          if(!nonParam) {warning("using nonparametric boot");nonParam<-TRUE}
          }
      if(is.list(o0)){
        est.psi00<-est.psi0<-o0$psi
        ss00<- o0$dev.no.gap
        if(!nonParam) fitted.ok<-fitted(o0)
        } else {
          if(!nonParam) stop("the first fit failed and I cannot extract fitted values for the semipar boot")
          if(random) {
            est.psi00<-est.psi0<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
            PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
            o0<-try(suppressWarnings(segConstr.glm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
            ss00<-o0$dev.no.gap
          } else {
          est.psi00<-est.psi0<-apply(PSI,2,mean)
          ss00<-opz$dev0
        }
        }

      n.intDev0<-nchar(strsplit(as.character(ss00),"\\.")[[1]][1])
      
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
      for(k in seq(n.boot)){
        #if(k==4) browser()
        ##se gli *ultimi* n.boot.rev valori di ss sono uguali, cambia i psi...
        n.boot.rev<- 3 #3 o 4?
        diff.selected.ss <- rev(diff(na.omit(all.selected.ss)))
        #if(length(na.omit(diff(all.selected.ss[1:n.boot.rev])))==(n.boot.rev-1) && all(round(diff(all.selected.ss[1:n.boot.rev]),6)==0)){
        if(length(diff.selected.ss)>=(n.boot.rev-1) && all(round(diff.selected.ss[1:(n.boot.rev-1)],6)==0)){
          qpsi<-sapply(1:ncol(Z),function(i)mean(est.psi0[i]>=Z[,i]))
          qpsi<-ifelse(abs(qpsi-.5)<.1, alpha, qpsi)
          alpha<-1-alpha
          est.psi0<-sapply(1:ncol(Z),function(i)quantile(Z[,i],probs=1-qpsi[i],names=FALSE))
        }
        
        #}
        PSI <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
        if(jt) Z<-apply(Z.orig,2,jitter)
        if(nonParam){
              id<-sample(n, size=size.boot, replace=TRUE)
              o.boot<-try(suppressWarnings(segConstr.glm.fit(y[id], XREG[id,,drop=FALSE], Z[id,,drop=FALSE], PSI[id,,drop=FALSE],
                w[id], offs[id], opz.boot)), silent=TRUE)
        } else {
              yy<-fitted.ok+sample(residuals(o0),size=n, replace=TRUE)
              o.boot<-try(suppressWarnings(segConstr.glm.fit(yy, XREG, Z.orig, PSI, weights, offs, opz.boot)), silent=TRUE)
        }
        if(is.list(o.boot)){
            all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi
        } else {
            est.psi.boot<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
        }
        #if(k==7) browser()
        ### se est.psi.boot non e' cambiato (e puoi vederlo da all.est.psi.boot), allora cambialo!
        
        
        PSI <- matrix(rep(est.psi.boot, rep(nrow(Z), length(est.psi.boot))), ncol = length(est.psi.boot))
        opz$h<-max(opz$h*.9, .2)
        opz$it.max<-opz$it.max+1
        o<-try(suppressWarnings(segConstr.glm.fit(y, XREG, Z.orig, PSI, w, offs, opz, return.all.sol=TRUE)), silent=TRUE)
        if(!is.list(o) && random){
                est.psi0<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
                PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
                o<-try(suppressWarnings(segConstr.glm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
                count.random<-count.random+1
        }
        #se il modello e' stato stimato controlla se la soluzione e' migliore..
        if(is.list(o)){
              if(!"coefficients"%in%names(o$obj)) o<-extract.psi(o)
              all.est.psi[k,]<-o$psi
              all.ss[k]<-o$dev.no.gap
              if(o$dev.no.gap<=ifelse(is.list(o0), o0$dev.no.gap, 10^12)) {o0<-o; k.psi.change<- k}
              est.psi0<-o0$psi
              all.selected.psi[k,] <- est.psi0
              all.selected.ss[k]<-o0$dev.no.gap #min(c(o$SumSquares.no.gap, o0$SumSquares.no.gap))
        }
            
        if(visualBoot) {
              flush.console()
              #      spp <- if (it < 10) " " else NULL
              #      cat(paste("iter = ", spp, it,
              #                "  dev = ",sprintf('%8.5f',L1), #formatC(L1,width=8, digits=5,format="f"), #era format="fg"
              #n.intDev0<-nchar(strsplit(as.character(dev.values[2]),"\\.")[[1]][1])
              cat(paste("boot sample = ", sprintf("%2.0f",k),
                        "  opt.dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), o0$dev.no.gap), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
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

#browser()
      
      # SS.ok<-min(all.selected.ss)
      # id.accept<- ((abs(all.ss-SS.ok)/SS.ok )<= 0.05)
      # psi.mean<-apply(all.est.psi[id.accept,,drop=FALSE], 2, mean)
      # est.psi0<-psi.mean
      # devi ristimare il modello con psi.mean
      # PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
      # o0<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)

      

      ris<-list(all.selected.psi=drop(all.selected.psi),all.selected.ss=all.selected.ss, all.psi=all.est.psi, all.ss=all.ss)

      if(is.null(o0$obj)){
          PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
          o0<-try(suppressWarnings(segConstr.glm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
      }
      if(!is.list(o0)) return(0)
      o0$boot.restart<-ris
      rm(.Random.seed, envir=globalenv())
      return(o0)
      }