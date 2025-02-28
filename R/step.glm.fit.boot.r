step.glm.fit.boot <- function(y, XREG, Z, PSI, w, offs, opz, n.boot=10, size.boot=NULL, jt=FALSE,
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
  #---------------------------------------------
  adj.psi <- function(psii, LIM) {
    pmin(pmax(LIM[1, ], psii), LIM[2, ])
  }
  #---
  extract.psi<-function(lista){
    #serve per estrarre il miglior psi..
    #dev.values<-lista[[1]][-1] #remove the 1st one referring to model without psi
    #psi.values<-lista[[2]][-1] #remove the 1st one (NA)
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
  #browser()
  if(is.null(opz$seed)){
    mY <- mean(as.numeric(y))
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
  visualBoot<-opz$display
  opz$display<-FALSE
  #opz.boot<-opz
  #opz.boot$pow=c(1,1) #c(1.1,1.2)
  opz1<-opz #opz1 viene usata solo quando diversi tentativi di stimare il modello falliscono...
  opz1$it.max <-0
  opz0 <- opz
  opz0$maxit.glm <- 2
  opz0$agg<-.2
  n<-length(y)
  alpha<-opz$alpha
  #limZ <- apply(Z, 2, quantile, names = FALSE, probs = alpha) #c(alpha, 1 - alpha))
  limZ <- if(is.null(opz$limZ)) apply(Z, 2, quantile, names=FALSE, probs=alpha) else opz$limZ
  rangeZ <- apply(Z, 2, range) #serve sempre
  o0 <-try(suppressWarnings(step.glm.fit(y, XREG, Z, PSI, w, offs, opz0, return.all.sol=FALSE)), silent=TRUE)
  #browser()
  if(!is.list(o0)) {
    o0<- suppressWarnings(step.glm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=TRUE))
    o0<-extract.psi(o0)
    ss00<-opz$dev0
    if(!nonParam) {warning("using nonparametric boot");nonParam<-TRUE}
  }
  if(is.list(o0)){
    est.psi00<-est.psi0<-o0$psi
    ss00<-o0$SumSquares.no.gap
    eta0 <- o0$eta0
    if(!nonParam) fitted.ok<-fitted(o0)
  } else {
    if(!nonParam) stop("the first fit failed and I cannot extract fitted values for the semipar boot")
    if(random) {
      est.psi00<-est.psi0<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
      PSI1 <- matrix(est.psi0, n, ncol = length(est.psi0), byrow=TRUE)
      o0<-try(suppressWarnings(step.glm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
      ss00<-o0$SumSquares.no.gap
      eta0 <- o0$eta0
    } else {
      est.psi00<-est.psi0<-apply(PSI,2,mean)
      ss00<-opz$dev0
      eta0 <- NULL
    }
  }
  
  n.intDev0<-nchar(strsplit(as.character(ss00),"\\.")[[1]][1])
  
  all.est.psi.boot<-all.selected.psi<-all.est.psi<-matrix(NA, nrow=n.boot, ncol=length(est.psi0))
  all.ss<-all.selected.ss<-rep(NA, n.boot)
  if(is.null(size.boot)) size.boot<-n
  
  Z.orig<-Z
  count.random<-0
  agg.values<-seq(.2,.05,l=n.boot)
  ###INIZIO BOOT
  alpha<-.1
  corr=1.2
  n.boot.rev<- 3 #3 o 4?
  
  for(k in seq(n.boot)){
    #if(k==7) browser()
    ##se gli *ultimi* n.boot.rev valori di ss sono uguali, cambia i psi...
    opz$eta0 <- eta0
    diff.selected.ss <- rev(diff(na.omit(all.selected.ss)))
    if(length(diff.selected.ss)>=(n.boot.rev-1) && all(round(diff.selected.ss[1:(n.boot.rev-1)],6)==0)){
      #browser()
      qpsi     <- sapply(1:ncol(Z),function(i)mean(est.psi0[i]>=Z[,i]))
      qpsi.cor <- sapply(1:ncol(Z),function(i)mean((corr*est.psi0[i])>=Z[,i]))
      qpsi <- ifelse(abs(qpsi-.5)<=.2, qpsi.cor, alpha)
      alpha<-1-alpha
      corr<-1/corr
      est.psi0 <- sapply(1:ncol(Z),function(i)quantile(Z[,i], probs=qpsi[i],names=FALSE))
      est.psi0 <- adj.psi(est.psi0, limZ)
      #est.psi0<- jitter(est.psi0, amount=min(diff(est.psi0))) 
    }
    ############################ 25/7/24 #####
    est.psi0 <- unlist(tapply(est.psi0, opz$id.psi.group, sort))
    ##########################################
    
    PSI <- matrix(est.psi0, n, ncol = length(est.psi0), byrow=TRUE)
    if(jt) Z<-apply(Z.orig,2,jitter)
    if(nonParam){
      id<-sample(n, size=size.boot, replace=TRUE)

      o.boot<-try(suppressWarnings(step.glm.fit(y[id], XREG[id,,drop=FALSE], Z[id,,drop=FALSE], PSI[id,,drop=FALSE],
                                            w[id], offs[id], opz)), silent=TRUE)
    } else {
      yy<-fitted.ok+sample(residuals(o0),size=n, replace=TRUE)
      o.boot<-try(suppressWarnings(step.glm.fit(yy, XREG, Z.orig, PSI, weights, offs, opz)), silent=TRUE)
    }
    if(is.list(o.boot)){
      all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi
    } else {
      est.psi.boot<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
      est.psi.boot<- unlist(tapply(est.psi.boot, opz$id.psi.group, sort))
    }
    PSI <- matrix(est.psi.boot, n, ncol = length(est.psi.boot), byrow=TRUE)
    #opz$h<-max(opz$h*.9, .2)
    opz$it.max<-opz$it.max+1
    opz$agg<-agg.values[k]
    #
    opz$Nboot <- k
    #
    o <-try(suppressWarnings(step.glm.fit(y, XREG, Z.orig, PSI, w, offs, opz, return.all.sol=TRUE)), silent=TRUE)
    #if(k==8) browser()
    
    if(!is.list(o) && random){
      est.psi0<-apply(limZ,2,function(r)runif(1,r[1],r[2]))
      PSI1 <- matrix(est.psi0, n, ncol = length(est.psi0), byrow=TRUE)
      o <-try(suppressWarnings(step.glm.fit(y, XREG, Z, PSI1, w, offs, opz1)), silent=TRUE)
      count.random<-count.random+1
    }
    #se il modello e' stato stimato controlla se la soluzione e' migliore..
    if(is.list(o)){
      if(!"coefficients"%in%names(o$obj)) o<-suppressWarnings(try(extract.psi(o), silent=TRUE))
      #if(class(o)!="try-error"){
      if(!inherits(o, "try-error")){
        #if(k==8) browser()
        all.est.psi[k,]<-o$psi[!is.na(o$psi)]
        all.ss[k]<- o$SumSquares.no.gap
        if(o$SumSquares.no.gap<=ifelse(is.list(o0), o0$SumSquares.no.gap, 10^12)) o0<-o
        est.psi0<-o0$psi
        all.selected.psi[k,] <- est.psi0
        all.selected.ss[k]<-L0<-o0$SumSquares.no.gap
        eta0 <- o0$eta0
      }
    }
    
    
    if (visualBoot) {
      flush.console()
      #      spp <- if (it < 10) " " else NULL
      #      cat(paste("iter = ", spp, it,
      #                "  dev = ",sprintf('%8.5f',L1), #formatC(L1,width=8, digits=5,format="f"), #era format="fg"
      #n.intDev0<-nchar(strsplit(as.character(dev.values[2]),"\\.")[[1]][1])
      unlpsi<- unlist(est.psi0)
      Lp<-length(unlpsi)
      
      cat(paste("boot sample = ", sprintf("%2.0f",k),
                "  opt.dev = ", sprintf(paste("%", n.intDev0+6, ".5f",sep=""), o0$SumSquares.no.gap), #formatC(L1,width=8, digits=5,format="f"), #era format="fg" 
                "  n.psi = ",formatC(Lp, digits=0,format="f"), 
                "  est.psi = ",paste(formatC(unlpsi[1:min(Lp,5)],digits=3,format="f"), collapse="  "), #sprintf('%.2f',x)
                sep=""), "\n")
    }
    #conta i valori ss uguali.. cosi puoi fermarti prima..
    asss<-na.omit(all.selected.ss)
    if(length(asss)>break.boot){
      if(all(rev(round(diff(asss),6))[1:(break.boot-1)]==0)) break
    }
  } #end n.boot
  
  all.selected.psi<-rbind(est.psi00,all.selected.psi)
  all.selected.ss<-c(ss00, all.selected.ss)
  
  #SS.ok<-min(all.selected.ss)
  #id.accept<- ((abs(all.ss-SS.ok)/SS.ok )<= 0.05)
  #psi.mean<-apply(all.est.psi[id.accept,,drop=FALSE], 2, mean)
  #est.psi0<-psi.mean
  # #devi ristimare il modello con psi.mean
  # PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
  # o0<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)
  
  
  
  ris<-list(all.selected.psi=drop(all.selected.psi),all.selected.ss=all.selected.ss, all.psi=all.est.psi, all.ss=all.ss)
  
  if(is.null(o0$obj)){
    PSI1 <- matrix(est.psi0, n, ncol = length(est.psi0),byrow=TRUE)
    o0 <- try(step.glm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)
    warning("The final fit can be unreliable (possibly mispecified stepmented relationship)", call.=FALSE, immediate.=TRUE)
  }
  if(!is.list(o0)) return(0)
  o0$boot.restart<-ris
  o0$seed <- seed
  
  #rm(.Random.seed, envir=globalenv())
  return(o0)
}
