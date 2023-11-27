selgmented <-function(olm, seg.Z, Kmax=2, type=c("score", "bic", "davies", "aic"), alpha=0.05, 
                  control=seg.control(), refit=FALSE, stop.if=6, return.fit=TRUE, bonferroni=FALSE, #improve.after.G=FALSE,
                  msg=TRUE, plot.ic=FALSE, th=NULL, G=1, check.dslope=TRUE){ #}, a=1){
  
  #vedere bene il fatto del th!!
  
  #check.dslope: if TRUE, a check on slope difference parameters is done to spot the "non-significant" ones. Then 
  #  the changepoints corresponding to such slopes are removed from the last fit. 
  #refit: if TRUE (and return.fit=TRUE), the last (selected) model is re-fitted using the 'control' controlling options (e.g. with n.boot>0)
  #stop.if. If type="bic" or "aic", the search of number of break can be stopped when the last 'stop.if' fits provide higher aic/bic value
  #th: When the distance between 2 estimated breakpoints is <=th, only one is retained. Default is th = drop(diff(rangeZ)/100)
  #mettere l'opzione di gdl=n.changepoint e poi (df*log(n))^alpha dove alpha=1.01 (vedi...)
  #sel1() si usa per G>1
  sel1 <- function(y, G, Kmax, type="bic", th=5, refit=FALSE, check.dslope=TRUE, msg=TRUE, bonferroni=FALSE){ #, a=1
    #BIC<-function(obj){
    #  n <-length(obj$residuals)
    #  r <- n*log(sum(obj$residuals^2)/n) + (n-obj$df.residual)*(log(n)^a) #- 1
    #  r
    #}
    BIC.f<-if(type=="bic") BIC else AIC
    #---
    drop.close <-function(all.psi, th){
      if(length(all.psi)==1) return(all.psi)
      all.psi <- c(1,sort(all.psi[!is.na(all.psi)]),n)
      id<- which(diff(all.psi)<=th)[1]
      while(!is.na(id) && length(id)>0){
        all.psi <- all.psi[-(id+1)]
        id<- which(diff(all.psi)<=th)[1]
      }
      #all.psi<- all.psi[-c(1, length(all.psi))]
      all.psi <- setdiff(all.psi,c(1, n))
      all.psi
    }
        #start..
    n<-length(y)
    x<-1:n
    n1<-ceiling(n/G)
    K1<-ceiling(Kmax/G)
    cutvalues <- c(seq(1,n,by=n1),n+1)
    id<-cut(x, cutvalues, right=FALSE, labels=FALSE)
    r<-vector("list",G)
    #browser()
    for(i in 1:G){
      yy<-y[id==i]
      xx<-x[id==i] 
      olm<-lm(yy~xx) #, data=d)
      .a <- capture.output(r[[i]] <-try(suppressWarnings(selgmented(olm, ~xx, type=type, Kmax=K1, 
                                                          refit=FALSE, msg=FALSE, G=1)), silent=TRUE)) 
      if(inherits(r[[i]],"try-error")) r[[i]]<- olm
      n.psi<- if(is.null(r[[i]]$psi) || any(is.na(r[[i]]$psi))) 0 else nrow(r[[i]]$psi)
      n.psi<- if(n.psi<10)  paste("", n.psi) else paste(n.psi)
      if(msg) {
        cat("\n##### subset ", paste(i, ": ...",sep=""))
        cat(" ",n.psi,"selected breakpoints \n")
      }
      
    }
    
    all.psi <-unlist(sapply(r, function(.x) .x$psi[,"Est."]))

    #browser()    
    
    psi.fromG <- drop.close(all.psi,th)
    
    psi.removed<- setdiff(all.psi, psi.fromG)
    if(length(psi.removed)>=1 && msg){
      cat(paste("\n", length(psi.removed), "breakpoint(s) removed for closeness (see argument 'th')\n"))
    }
    
    all.psi <- psi.fromG


    olm <-lm(y~x)
    
    newpsi <- cutvalues[-c(1, length(cutvalues))]
    if(msg){
      cat(paste("\n Assessing the", length(newpsi), "cutoff(s) as breakpoint(s)..\n"))
    }
    
    all.psi<-sort(c(all.psi, newpsi)) #tutti i psi
    
    #elimina quelli vicini..
    psi.withCut <- all.psi <-drop.close(all.psi,th)
    
    
    ######## ORA stima modelli con un numero decrescente di psi...
    tvalueU <- psi0 <- all.psi
    list.psi <- list.fit<-NULL

    #browser()
    #se nella riduzione del numero di  psi il modello non viene stimato, allora piuttosto che passargli i valori di psi
    # da cui ha arbitrariamente eliminato il primo, passagli il numero di psi..
    fit.ok<-TRUE
    while(length(tvalueU)>=1){
      if(fit.ok){
        .a <- capture.output(os0 <- try(suppressWarnings(
          segmented(olm, ~x, psi=all.psi, control=seg.control(n.boot=0, alpha=.01))), silent=TRUE))
      } else {
        .a <- capture.output(os0 <- try(suppressWarnings(
          segmented(olm, ~x, npsi=length(all.psi), control=seg.control(n.boot=0, alpha=.01))), silent=TRUE))
      }
      if(inherits(os0, "segmented")) {
        fit.ok<-TRUE
        list.fit[[length(list.fit)+1]] <- os0
        list.psi[[length(list.psi)+1]]<- os0$psi[,"Est."]
        if(length(os0$psi[,"Est."])==1) break
        tvalueU<- abs(summary(os0)$coefficients[os0$nameUV$U,3])
        idU <- which.min(abs(tvalueU))
        all.psi <- os0$psi[-idU,"Est."]
      } else {
        fit.ok<-FALSE
        list.psi[[length(list.psi)+1]]<-list.fit[[length(list.fit)+1]] <- NA
        tvalueU<-all.psi<-all.psi[-1]
      }
    
    }
    #browser()
    bicV <-  sapply(list.fit, function(.x) if(inherits(.x, "segmented")) BIC(.x) else NA)
    if(length(bicV)!= length(psi.withCut)) stop("Errore nella dim!!")
    list.fit[[length(list.fit)+1]] <- olm
    bicV<- c(bicV, BIC.f(olm))
    r <- r0 <- list.fit[[which.min(bicV)]]
    npsi.ok <- ((length(psi.withCut)):0)[which.min(bicV)]
    
    #DUBBIO: il controllo ed eliminazione dei psi lo facciamo all'interno del "while(length(tvalueU)>1)"?
    #   ed inoltre al modello selezionato un ultimo fit con boot bisognerebbe darlo...
    #browser()
    
    if(inherits(r, "segmented")){
      #ATTENZIONE QUA CI VORREBBE UN ALTRO CONTROLLO SULLA VICINANZA DEI psi..
      all.psi <- sort(r$psi[,"Est."])
      #all.psi <- drop.close(all.psi,th)
      #r$psi<-all.psi
      #browser()
      cont1 <- length(all.psi)>0 && (length(all.psi)!=length(drop.close(all.psi,th))) 
        #length(all.psi)!=length(drop.close(all.psi,th))
      while(cont1){
        start.psi<-drop.close(all.psi,th)
        r0$call$psi<- start.psi 
        .a <- capture.output(r <- suppressWarnings(try(update(r0))), type="message")
        all.psi<- if(inherits(r, "segmented")) r$psi[,2] else start.psi[-1]
        cont1<- !(inherits(r, "segmented") && (length(all.psi)==length(drop.close(all.psi,th))))
        cont1 <- cont1 &&  length(all.psi)>0
      }
      
      if(inherits(r, "segmented")){
      if(check.dslope){
        all.psi<- r$psi[,"Est."]
        #rm.id <- which(abs(slope(r)[[1]][,3]) <= qnorm(1-alpha/2))
        soglia <- if(!bonferroni) qnorm(1-alpha/2) else qnorm(1-alpha/(2* length(r$nameUV$U)) ) 
        rm.id <- which(abs(summary.lm(r)$coefficients[r$nameUV$U, 3]) <= soglia) #era "t value" invece che 3 ma non funzionava con glm..
          
        if(length(rm.id)>0){
          all.psi <- all.psi[-rm.id]
          if(length(all.psi)>0){
            if(refit) {
              control$alpha <- .05
              r$call$control<-control
            }
            r$call$psi=all.psi
            #.a <- capture.output(r<-suppressWarnings(try(r<-update(r), silent=TRUE)))
            .a <- capture.output(r <- try(suppressWarnings(update(r))))
            if(!inherits(r,"segmented")){
              #facciamo un altro tentativo..
              r0$call$psi=all.psi
              control$alpha <- .005
              r0$call$control<-control
              .a <- capture.output(r<-suppressWarnings(try(r<-update(r0), silent=TRUE)))
              #r<-update(r0)
            }
            n.psi.ok <- if(is.null(r$psi)) 0 else nrow(r$psi)
          } else {
            r <- olm
            n.psi.ok <- 0  
          }
        } else {
          n.psi.ok <- nrow(r$psi)
          }
      } else {
        n.psi.ok <- nrow(r$psi)
      }
    } else {n.psi.ok <- 0}    
} else {
      n.psi.ok <- 0
    }
    if(msg) cat("\n####### Overall: ...  ", n.psi.ok, "selected breakpoint(s) \n\n")
    #cat(" ##### Overall ", nrow(r$psi),"selected breakpoints \n")
    if(!is.list(r)) r<- olm
    r$cutvalues <- cutvalues #including the extremes
    r$selection.psi <-list(bic.values=bicV, npsi=n.psi.ok)
    r
    #r<-segmented(olm, ~x, psi=r$psi[,"Est."], control=seg.control(n.boot=10,alpha=.01)) #fix.npsi=FALSE
    #r
  } #fine sel1()
  #=====================================================================
  BIC<-function(obj, a=1){
    n <-length(obj$residuals)
    r <- n*log(sum(obj$residuals^2)/n) + (n-obj$df.residual)*(log(n)^a) #- 1
    r
  }
  #=====================================================================
  type<-match.arg(type)
  if(!type%in%c("bic","aic") && Kmax!=2) stop("Kmax>2 is not (yet?) allowed with hypothesis testing procedures", call.=FALSE)
  if(!type%in%c("bic","aic") && G>1) stop("G>1 is allowed only with type='aic' or 'bic' ", call.=FALSE)
  #=====================================================================
  if(G==1){  
    build.mf<-function(o, data=NULL){
      #returns the dataframe including the possibly untransformed variables,
      #including weight and offset
      fo<-formula(o)
      if(!is.null(o$weights))
        fo<-update.formula(fo,paste("~.+",all.vars(o$call$weights), sep="")) 
      if(!is.null(o$call$offset))
        fo<-update.formula(fo,paste("~.+",all.vars(o$call$offset), sep="")) 
      if(!is.null(o$call$subset))
        fo<-update.formula(fo,paste("~.+",all.vars(o$call$subset), sep=""))
      #o$call$formula<-fo
      if(is.null(o$call$data)) {
        R<-get_all_vars(fo)
      } else { 
        R<-get_all_vars(fo, data=eval(o$call$data))
      }
      R
    }
    
    #browser()
    
    if(is.numeric(olm)){
      y<-olm
      Z<-x<- 1:length(y)
      olm <- lm(y~x)
    } else {
      if(missing(seg.Z)){
        nomeX <- all.vars(formula(olm))[2]
        if(length(nomeX)>1 || any(is.na(nomeX))) stop("I cannot determine the segmented variable")
        seg.Z<- as.formula(paste("~", nomeX ))
        Z <- olm$model[[nomeX]]
      } else {
        if(length(all.vars(seg.Z))>1) stop("Multiple variables are not allowed in seg.Z")
      }
    }
    if(type%in%c("bic","aic")){
      control1<-control
      control1$n.boot = 0
      control1$tol <- .001 #default e' .00001
      control1$alpha<-.01
      ICname<- if(type=="bic") "BIC" else "AIC"
      BIC.f<-if(type=="bic") BIC else AIC
      bicM0 <- BIC.f(olm)
      
      Kmax <- min(floor((olm$df.residual-1)/2), Kmax)
      
      npsi<-1:Kmax
      startpsi<-ris<-vector("list", length(npsi))  
      conv<-bic.values<- rep(NA, length(npsi))
      
      if(!is.null(olm$call$data)) assign(paste(olm$call$data), eval(olm$call$data, envir=parent.frame() ))
      
      #fit with 1 breakpoint
      .a<-capture.output(os<- suppressWarnings(try(segmented(olm, seg.Z, npsi=1, control=control1), silent=TRUE)))
      ris[[1]] <- os #<- suppressWarnings(try(segmented(olm, seg.Z, npsi=1, control=control1), silent=TRUE))
      #if fails try boot restating
      if(inherits(os, "try-error")) {
        .a <- capture.output(os<- suppressWarnings(try(segmented(olm, seg.Z, npsi=1, control=control), silent=TRUE)))
        ris[[1]]<- os
      }
      
      if(inherits(os, "segmented")){
        if(is.null(th)) th <- drop(diff(os$rangeZ)/100)
        bic.values[1]<- BIC.f(os)
        Z<- os$model[,os$nameUV$Z]
        m1 <-min(Z)
        m2 <-max(Z)
        estpsi <- os$psi[,"Est."]
        M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
        psi0 <- sum(M[which.max(apply(M,1,diff)),])/2
        startpsi[[1]] <- sort(c(estpsi, psi0))
        conv[1] <- 1 
      } else {
        m1 <-min(Z)
        m2 <-max(Z)
        estpsi <- (m1+m2)/2
        M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
        #psi0 <- sum(M[which.max(apply(M,1,diff)),])/2
        startpsi[[1]] <- estpsi
        bic.values[1]<- BIC.f(olm)+1
        conv[1] <- 0
      }
      
      i=1 #ponilo =1 
      if(Kmax>=2){
        if(msg) {
          flush.console()
          cat(paste("No. of breakpoints: "))
        }
        
        for(i in 2:Kmax){
          #if(i==14) browser()
          .a <- capture.output(os<-suppressWarnings(try(segmented(olm, seg.Z, psi=startpsi[[i-1]], control=control1), silent=TRUE)))
          ris[[i]]<- os
          if(msg) {
            flush.console()
            cat(paste(i,".. "))
          }
          
          if(inherits(os, "segmented")) {
            conv[i]<-1
            estpsi <- ris[[i]]$psi[,"Est."]
            id<- which(diff(estpsi)<=th)
            if(length(id)>0){
              estpsi <- estpsi[-id]
            }
            M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
            diffpsi <- apply(M,1,diff)
            if(length(id)>0){
              psi0<-NULL
              for(j in 1:(length(id))) {
                psi0[length(psi0)+1] <- sum(M[which(diffpsi==rev(sort(diffpsi))[j]),])/2
              }
            } else {
              id.max.diffpsi <- which.max(diffpsi)
              psi0 <- sum(M[id.max.diffpsi,])/2
            }
            startpsi[[i]] <- sort(c(estpsi, psi0))
            bic.values[i]<- BIC.f(ris[[i]]) #-2*logLik(ris[[i]]))+ edf*log(n)*Cn
          } else {
          #se non e' arrivato a convergenza comunque aggiungi un breakpoint..
          #   SENSATO???x
            control1$alpha<-.005
            conv[i]<-0
            bic.values[i]<- bic.values[i-1]+1
            M<-matrix(c(m1 ,rep(startpsi[[i-1]], each=2), m2), 
                      ncol=2, nrow=length(startpsi[[i-1]])+1, byrow=TRUE)
            diffpsi <- apply(M,1,diff)
            #psi0 <- sum(M[which.max(diffpsi[-id.max.diffpsi]),])/2
            psi0 <- sum(M[which.max(diffpsi),])/2
            startpsi[[i]] <- sort(c(psi0, M[-c(which.min(diffpsi), nrow(M)),2]))
            #startpsi[[i]] <- sort(c(estpsi, psi0))
          }
          #un controllo su bic values.. fermarsi SE gli ultimi K sono NA oppure sono crescenti!!
          if(i>=stop.if && all(rev(na.omit(diff(c(bicM0,bic.values))))[1:stop.if]>=0)) break
        } #end in for(2 in Kmax)
      
        if(msg) cat("\n")  
        }
      #browser()
      
      npsiVeri <-sapply(ris, function(.x) if(is.list(.x))nrow(.x$psi) else NA )[1:i]
      #npsiVeri <-c(1,unlist(sapply(startpsi, length))[1:(i-1)])
      
      bic.valuesOrig <- bic.values 
      bic.values<- bic.values[1:i] #bic.values[!is.na(bic.values)]
      bic.values <- c(bicM0, bic.values)
      names(bic.values)<-c("0", npsiVeri)
      n.psi.ok<- c(0,npsiVeri)[which.min(bic.values)]

      if(n.psi.ok==0){
        m<-matrix(NA,1,1, dimnames=list(NULL, "Est."))
        olm$selection.psi<- list(bic.values=bic.values, npsi=n.psi.ok)
        olm$psi<-m
        if(msg){
          if(any(is.na(bic.valuesOrig))) cat(paste("(search truncated at ", i, " breakpoints due to increasing values of ",  ICname ,") \n", sep=""))
          cat(paste("\n",ICname, " to detect no. of breakpoints:\n",sep=""))
          print(bic.values)
          cat(paste("\nNo. of selected breakpoints: ", n.psi.ok, " \n"))
        }
        return(olm)
      }
      
      if(i==n.psi.ok && msg) warning(paste("The best",ICname, "value at the boundary. Increase 'Kmax'?"), call.=FALSE, immediate. = TRUE)
      id.best<- which.min(bic.values[-1])
      
      if(!return.fit) {
        r <- list(bic.values=bic.values, npsi=n.psi.ok)
        return(r) 
      }
      #browser()
      r<- r0 <- ris[[id.best]]
      
      f<-function(x, soglia){
        #restituisce l'indice del vettore x t.c. il valore e' il piu' piccolo 
        #tra quelli che sono minori della soglia  
        id <- (x<=soglia)
        xx <- x[id]
        ind <- (1:length(x))[id]
        id.ok <- which.min(xx)
        ind[id.ok]
      }
      
      rm.after.check <- 0
      if(check.dslope){
        all.psi<- r$psi[,"Est."]
        soglia <- if(!bonferroni) qnorm(1-alpha/2) else qnorm(1-alpha/(2* length(r$nameUV$U))) 
        tU <- abs(summary(r)$coefficients[r$nameUV$U, 3])
        #if(length(tU[tU<=soglia])==length(tU)) #anche se tutti i t<= soglia fai comunque la procedura, perche' 
        #riducendo i psi, i tU potrebbero cambiare
        #rm.id <- which.min(tU[tU<=soglia])
        rm.id <- f(tU, soglia)
        while(length(rm.id)>0){
          rm.after.check <- rm.after.check+1
          all.psi <- all.psi[-rm.id]
          if(length(all.psi)<=0) break
          r0$call$psi=quote(all.psi)
          .a <- capture.output(r <- suppressWarnings(try(update(r0), silent=TRUE)))
          if(!inherits(r,"segmented")){ #facciamo un altro tentativo..
            #r0$call$psi=all.psi
            control$alpha <- .005
            r0$call$control<-control
            .a <- capture.output(r<-suppressWarnings(try(r<-update(r0), silent=TRUE)))
          }
          if(inherits(r,"segmented")){
            tU <- abs(summary(r)$coefficients[r$nameUV$U, 3])
            #rm.id <- which.min(tU[tU<=soglia])
            rm.id <- f(tU, soglia)
            all.psi<-r$psi[,"Est."]
          } else {
            rm.id<-1
          }
        }
        if(length(all.psi)<=0 || is.null(r$psi) ) {
          n.psi.ok<-0
          r<- olm
        } else {
          n.psi.ok <-  nrow(r$psi)
        }
      } else { #end if(check.dslope)
        #ATTENZIONE: SE NON HA SELEZIONATO BREAKPOINTS???
    #browser()
        if(refit){
          r$call$psi<- r$psi[,"Est."]#startpsi[[n.psi.ok-1]]
          control$alpha <- .005
          r$call$control<-control
          r <- update(r)
        }
        #browser()
        n.psi.ok<-length(r$psi[,2]) #in realta' gia' c'e' "n.psi.ok"
        
      }
      if(plot.ic) {
        plot(c(0,npsiVeri), bic.values, xlab=" No. of breakpoints", ylab=ICname, type="o")
        points(n.psi.ok, min(bic.values), pch=19)
      }
      #browser()
      if(msg){
        if(any(is.na(bic.valuesOrig))) cat(paste("(search truncated at ", i, " breakpoints due to increasing values of ",  ICname ,") \n", sep=""))
        cat(paste("\n",ICname, " to detect no. of breakpoints:\n",sep=""))
        print(bic.values)
        add.msg <- if(rm.after.check==0) " \n" else paste(" (", rm.after.check, " breakpoint(s) removed due to small slope diff)\n", sep="")
        cat(paste("\nNo. of selected breakpoints:", n.psi.ok, add.msg))
      }
      r$selection.psi <- list(bic.values=bic.values, npsi=n.psi.ok)
      return(r)
    } #end aic/bic..
    alpha.adj<-alpha/Kmax
    p1<- if(type=="score")  pscore.test(olm, seg.Z, n.break=2)$p.value else davies.test(olm)$p.value
    p1.label<-"p-value '0 vs 2' "
    if(p1>alpha.adj){
      p2.label<-"p-value '0 vs 1' "
      p2<- if(type=="score") pscore.test(olm, seg.Z, n.break=1)$p.value else p1 #davies.test(olm)$p.value
      if(!bonferroni) alpha.adj<- alpha
      if(p2>alpha.adj) {
        out<-olm
      } else {
        out<-segmented(olm, seg.Z, npsi=1, control=control)
      }
    } else {
      p2.label<-"p-value '1 vs 2' "
      #################
      #browser()
      #MF<-build.mf(olm)
      #olm<-update(olm, data=MF)
      #olm$call$data<-quote(MF)
      
      #olm<-update(olm, data=model.frame(olm)) #questo e' necessario per far funzionare davies.test() sotto..
      ################
      
      if(type=="score") {
        o1<-segmented(olm, seg.Z, npsi=1, control=control)
        p2<-pscore.test(o1, seg.Z, more.break=TRUE)$p.value
      } else {
        #KK<-new.env()
        #olm1<-update(olm, data=model.frame(o1))
        #o1<-  update(o1, obj=olm1)
        MF<-build.mf(olm)
        olm<-update(olm, data=MF)
        # olm$call$data<-quote(MF)
        #olm<-update(olm, data=model.frame(olm)) #questo e' necessario per far funzionare davies.test() sotto..
        
        o1 <- segmented(olm, seg.Z, npsi = 1, control = control)
        p2<-  davies.test(o1, seg.Z)$p.value
      }
      if(!bonferroni) alpha.adj<-alpha 
      if(p2>alpha.adj) {
        o1<-segmented(olm, seg.Z, npsi=1, control=control)
        #cat("One breakpoint detected\n")
        out<-o1
      } else {
        o2<-segmented(olm, seg.Z, npsi=2, control=control)
        #cat("Two breakpoint detected\n")
        out<-o2
      }
    }
    n.psi.ok<-length(out$psi[,2])
    x2<- -2*sum(log(c(p1,p2)))
    p<-1-pchisq(x2, df=2*2)
    r<-list(pvalues=c(p1=p1, p2=p2, p=p), npsi=n.psi.ok)
    attr(r, "label")<- p2.label
    if(!return.fit) {
      return(r)
    }
    if(msg){
      cat("Hypothesis testing to detect no. of breakpoints\n")
      type <- chartr(strsplit(type,"")[[1]][1], toupper(strsplit(type,"")[[1]][1]), type) #serve per render maiuscola la prima lettera..
      cat(paste("statistic:", type,"  level:", alpha, "  Bonferroni correction:", bonferroni, "\n"))
      cat(paste(p1.label, "= ", format.pval(p1,4), "   ", p2.label, "= ", format.pval(p2,4) ,
                " \nOverall p-value = ", format.pval(p,4),"\n",sep=""))
      cat(paste("No. of selected breakpoints: ", n.psi.ok, "\n"))
    }
    out$selection.psi<-r
    return(out)
  } else { #se G>1
    #browser()
    if(!is.numeric(olm)) stop("If G>1, 'olm' should be a numeric vector")
    if(is.null(th)) th <- max(round(length(olm)/100), 5)
    r <- sel1(y=olm, G=G, Kmax=Kmax, type=type, th=th, refit=FALSE, check.dslope = check.dslope, 
              msg=msg, bonferroni=bonferroni)
    r
  }
}

