selgmented <-function(olm, seg.Z, Kmax=2, type=c("score", "bic", "davies", "aic"), alpha=0.05, 
                  control=seg.control(), refit=FALSE, stop.if=5, return.fit=TRUE, bonferroni=FALSE, #improve.after.G=FALSE,
                  msg=TRUE, plot.ic=FALSE, th=NULL, G=1, check.dslope=TRUE){ #}, a=1){
  
  #vedere bene il fatto del th!!
  
  #ruolo di alpha??? e se la ricerca si vuole restringere a un intervallo? perndere seg.control()$alpha? 
  
  #check.dslope: if TRUE, a check on slope difference parameters is done to spot the "non-significant" ones. Then 
  #  the changepoints corresponding to such slopes are removed from the last fit. 
  #refit: if TRUE (and return.fit=TRUE), the last (selected) model is re-fitted using the 'control' controlling options (e.g. with n.boot>0)
  #stop.if. If type="bic" or "aic", the search of number of break can be stopped when the last 'stop.if' fits provide higher aic/bic value
  #th: When the distance between 2 estimated breakpoints is <=th, only one is retained. Default is th = drop(diff(rangeZ)/100)
  #mettere l'opzione di gdl=n.changepoint e poi (df*log(n))^alpha dove alpha=1.01 (vedi...)
  #sel1() si usa per G>1
  #===
  if(stop.if<=0) stop("'stop.if' should be an integer (at least 4, probably)")
  stop.if<-ceiling(stop.if)
  f<-function(x, soglia){
    #restituisce l'indice del vettore x t.c. il valore e' il piu' piccolo 
    #tra quelli che sono minori della soglia  
    id <- (x<=soglia)
    xx <- x[id]
    ind <- (1:length(x))[id]
    id.ok <- which.min(xx)
    ind[id.ok]
  }
  #===
  sel1 <- function(y, x, G, Kmax, type="bic",th=th,refit=FALSE,check.dslope=TRUE,
                   msg=TRUE, bonferroni=FALSE, olm0, control){ #, a=1
    #BIC<-function(obj){
    #  n <-length(obj$residuals)
    #  r <- n*log(sum(obj$residuals^2)/n) + (n-obj$df.residual)*(log(n)^a) #- 1
    #  r
    #}
    ICname<- if(type=="bic") "BIC" else "AIC"
    
    control1<-control
    control1$n.boot = 0
    control1$tol <- .001 #
    
    
    BIC.f<-if(type=="bic") BIC else AIC
    #---
    drop.close <-function(all.psi, th){
      if(length(all.psi)==1) return(all.psi)
      all.psi <- c(m1, sort(all.psi[!is.na(all.psi)]), m2)
      id<- which(diff(all.psi)<=th)[1]
      while(!is.na(id) && length(id)>0){
        all.psi <- all.psi[-(id+1)]
        id<- which(diff(all.psi)<=th)[1]
      }
      #all.psi<- all.psi[-c(1, length(all.psi))]
      all.psi <- setdiff(all.psi,c(m1, m2))
      all.psi
    }
        #start..
    n<-length(y)
    K1<-ceiling(Kmax/G)
    
    #browser()
    #x<-1:n
    #n1<-ceiling(n/G)
    #cutvalues <- c(seq(1,n,by=n1),n+1)
    
    cutvalues <- c(min(x), cumsum(rep(sum(range(x))/G, G)))

    id<-cut(x, cutvalues, right=FALSE, labels=FALSE)
    r<-vector("list",G)
    
    #browser()
    
    for(i in 1:G){
      #if(i==4) browser()
      yy<-y[id==i]
      xx<-x[id==i] 
      olm<-lm(yy~xx) #, data=d)
      .a <- capture.output(r[[i]] <-try(suppressWarnings(selgmented(olm, ~xx, type=type, Kmax=K1, 
                                                    refit=FALSE, msg=FALSE, G=1, control=control, check.dslope=FALSE)
                                                    ), silent=TRUE)) 
      if(inherits(r[[i]],"try-error")) r[[i]]<- olm
      n.psi<- if(is.null(r[[i]]$psi) || any(is.na(r[[i]]$psi[,"Est."]))) 0 else nrow(r[[i]]$psi)
      n.psi<- if(n.psi<10)  paste("", n.psi) else paste(n.psi)
      if(msg) {
        cat("\n##### subset ", paste(i, ": ...",sep=""))
        cat(" ",n.psi,"selected breakpoints \n")
      }
      
    }
    
    #browser()
    
    all.psi <-unlist(sapply(r, function(.x) .x$psi[,"Est."]))
    psi.fromG <- drop.close(all.psi,th)
    psi.removed<- setdiff(all.psi, psi.fromG)
    psi.removed<- psi.removed[!is.na(psi.removed)]
    if(length(psi.removed)>=1 && msg){
      cat(paste("\n", length(psi.removed), "breakpoint(s) removed for closeness (see argument 'th')\n"))
    }
    
    all.psi <- psi.fromG

    olm <- olm0 #lm(y~x)
    newpsi <- cutvalues[-c(1, length(cutvalues))]
    if(msg){
      cat(paste("\n => Assessing the", length(newpsi), "cutoff(s) as breakpoint(s). Computing the", ICname, "values.. \n"))
    }
    
    all.psi<-sort(c(all.psi, newpsi)) #tutti i psi
    
    #elimina quelli vicini..
    psi.withCut <- all.psi <-drop.close(all.psi,th)

    
    #browser()
    
    ######## ORA stima modelli con un numero decrescente di psi...
    tvalueU <- psi0 <- all.psi
    list.psi <- list.fit<-NULL

    #browser()
    #se nella riduzione del numero di  psi il modello non viene stimato, allora piuttosto che passargli i valori di psi
    # da cui ha arbitrariamente eliminato il primo, passagli il numero di psi..
    #fit.ok<-TRUE
    idbreak<-FALSE
    conv <- bicVa <- NULL
    #browser()
    if(msg) cat("    no. breakpoints: ", length(all.psi))
    control2<-control1
    control2$n.boot=6
    while(length(tvalueU)>=1){
      .a <- capture.output(os0 <- try(suppressWarnings(
          segmented(olm, ~x, psi=all.psi, control=control2)), silent=TRUE))
      if(!inherits(os0, "segmented")) {
        .a <- capture.output(os0 <- try(suppressWarnings(
          segmented(olm, ~x, npsi=length(all.psi), control=control1)), silent=TRUE))
      }
      
      if(inherits(os0, "segmented")) {
        conv[length(conv)+1]<- 1
        bicVa[length(bicVa)+1] <- BIC.f(os0)
        list.fit[[length(list.fit)+1]] <- os0
        list.psi[[length(list.psi)+1]]<- os0$psi[,"Est."]
        if(length(os0$psi[,"Est."])==1) break
        tvalueU<- abs(summary(os0)$coefficients[os0$nameUV$U,3])
        idU <- which.min(abs(tvalueU))
        #man mano che rimuovi i psi, se il t della diffSlope del psi che stai rimuovendo e' > soglia, fermati!
        ## allora non continuare a toglierli..
        #===>>>NOOOOO!! puo' accadere poi che toglie qualche psi, si ferma perche' tutti i t sono grandi, ma il bic e' piu' basso
        #in un modello precedente.. Quindi alla fine il criterio non e' ne bic e neanche tutti i t significativi.
        #Allora, coerentemente con il caso di G=1, la "verifica delle slope nonsignif, va fatta sul modello 
        #selezionato con il BIC!! Quindi commentiamo le righe di sotto (il idbreak si potrebbe pure eliminare)
        #e il controllo checkdslope lo facciamo dopo sul modello selezionato dal bic.
        #browser()
        #if(check.dslope){
        #  soglia <- if(!bonferroni) qnorm(1-alpha/2) else qnorm(1-alpha/(2* length(os0$nameUV$U)))
        #  if(abs(tvalueU[idU])>soglia) {idbreak=TRUE; break}
        #}
        all.psi <- os0$psi[-idU,"Est."]
      } else {
        conv[length(conv)+1]<- 0
        list.psi[[length(list.psi)+1]]<-list.fit[[length(list.fit)+1]] <- NA
        tvalueU<-all.psi<-all.psi[-1]
        bicVa[length(bicVa)+1] <- if(is.numeric(bicVa[length(bicVa)])) {bicVa[length(bicVa)]+1} else {1e4}
      }
      if(msg) cat(" ..", length(all.psi))
      #browser()
      if(length(bicVa)>stop.if && all( na.omit(rev(diff(bicVa))[1:stop.if])>0) && (sum(na.omit(rev(conv)[1:stop.if]))>0)) {idbreak<-TRUE; break}
      #na.omit() sta per eliminare gli NA che si creano se diff(bicVa) ha dimensione < stop.if
      #la sum(na.omit(rev(conv)[1:stop.if]))>0 indica che il controllo sui valori del bic va fatto solo 
      #se gli ultimi fit non sono tutti insuccessi. 
    }
    if(msg) cat(" .. 0\n")
    bicV <-  sapply(list.fit, function(.x) if(inherits(.x, "segmented")) BIC.f(.x) else NA)
    
    if(length(bicV)!=length(bicVa)) stop("Errore inatteso 1")
    
    bicV<- bicVa
    
    #browser()
    
    if(idbreak){
      if(msg) cat(paste(" =>", length(psi.withCut)-length(bicVa), "unevaluated model(s) due to", stop.if, "increasing A/BIC value(s)..\n"))
      #se si ferma prima, significa che sono stati valutati length(bicV) modelli con numero di psi da 
      #"length(psi.withCut)" fino a "length(psi.withCut)-length(bicV)+1" 
      #length(psi.withCut)-length(bicV)+1
      r <- r0 <- list.fit[[which.min(bicV)]]
      npsi.ok <-  (length(psi.withCut):(length(psi.withCut)-length(bicV)+1))[which.min(bicV)]
      nameBIC<-paste((length(psi.withCut):(length(psi.withCut)-length(bicV)+1)))
    } else {
      if(length(bicV)!= length(psi.withCut)) stop("Errore nella dim!!")
      list.fit[[length(list.fit)+1]] <- olm
      bicV<- c(bicV, BIC.f(olm))
      r <- r0 <- list.fit[[which.min(bicV)]]
      npsi.ok <- ((length(psi.withCut)):0)[which.min(bicV)]
      if(npsi.ok!=length(r$psi[,"Est."])) stop("Unexpected error..")
      nameBIC<-paste((length(psi.withCut)):0)
    }
    names(bicV)<-nameBIC
    
    #DUBBIO: il controllo ed eliminazione dei psi lo facciamo all'interno del "while(length(tvalueU)>1)"?
    #   ed inoltre al modello selezionato un ultimo fit con boot bisognerebbe darlo...
    #browser()
    
    if(inherits(r, "segmented")){
      #ALTRO CONTROLLO SULLA VICINANZA DEI psi..
      all.psi <- sort(r$psi[,"Est."])
      #all.psi <- drop.close(all.psi,th)
      #r$psi<-all.psi
      #browser()
      cont1 <- length(all.psi)>0 && (length(all.psi)!=length(drop.close(all.psi,th))) 
      while(cont1){
        start.psi<-drop.close(all.psi,th)
        r0$call$psi<- start.psi 
        .a <- capture.output(r <- suppressWarnings(try(update(r0))), type="message")
        all.psi<- if(inherits(r, "segmented")) r$psi[,"Est."] else start.psi[-1]
        cont1<- !(inherits(r, "segmented") && (length(all.psi)==length(drop.close(all.psi,th))))
        cont1 <- cont1 &&  length(all.psi)>0
      }
      if(inherits(r, "segmented")){
        # all.psi<- r$psi[,"Est."]
        # soglia <- if(!bonferroni) qnorm(1-alpha/2) else qnorm(1-alpha/(2* length(r$nameUV$U)) ) 
        # rm.id <- which(abs(summary(r)$coefficients[r$nameUV$U, 3]) <= soglia)
        # #===============
        # if(check.dslope && length(rm.id)>0){
        #   #     se devi controllare le slopeDiff
        #   #     #ALTRO CONTROLLO SULLA SIGNIF DELLE SLOPE-DIFF 
        #   #Non puoi eliminarli tutti assieme, perche' una volta che il piu' piccolo e stato eliminato 
        #gli altri possono cambiare..
          #browser()
        
        rm.after.check <- 0
        all.psi<- r$psi[,"Est."]
        #browser()
        if(check.dslope){
          soglia <- if(!bonferroni) qnorm(1-alpha/2) else qnorm(1-alpha/(2* length(r$nameUV$U))) 
          tU <- abs(summary(r)$coefficients[r$nameUV$U, 3])
          rm.id <- f(tU, soglia)
          while(length(rm.id)>0){
            rm.after.check <- rm.after.check+1
            all.psi <- all.psi[-rm.id]
            if(length(all.psi)<=0) break
            r0$call$psi=quote(all.psi)
            .a <- capture.output(r <- suppressWarnings(try(update(r0), silent=TRUE))) #senza boot
            if(!inherits(r,"segmented")){ #facciamo un altro tentativo..
              #r0$call$psi=all.psi
              #control$alpha <- .005
              r0$call$control<-control #con boot
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
            if(msg) warning("All found psi had a non-significant slope diff and have been removed", call.=TRUE, immediate.=TRUE)
            n.psi.ok<-0
            r<- olm
          } else {
            if(msg) cat(paste(" => ", rm.after.check, " breakpoint(s) removed due to 'small' slope difference\n", sep=""))
            if(refit){
              all.psi<- r$psi[,"Est."]
              r$call$psi=quote(all.psi)
              r$call$control<- quote(control)
              .a <- capture.output(r <- try(suppressWarnings(update(r0))))
              n.psi.ok <- if(is.null(r$psi)) 0 else nrow(r$psi)
            } else {
              all.psi<- r$psi[,"Est."]
              n.psi.ok <- nrow(r$psi)
            }
          }
        } else {#end if(check.dslope)
          
          if(refit){
            all.psi<- r$psi[,"Est."]
            r$call$psi=quote(all.psi)
            r$call$control<- quote(control)
            .a <- capture.output(r <- try(suppressWarnings(update(r))))
            n.psi.ok <- if(is.null(r$psi)) 0 else nrow(r$psi)
          } else {
            all.psi<- r$psi[,"Est."]
            n.psi.ok <- nrow(r$psi)
          }
        }
      } else {
        #stop("Errore inatteso 1")
       r<-olm
       n.psi.ok<-0
      }
    } else {
      #stop("Errore inatteso 2")
      r<-olm
      n.psi.ok<-0
    }

    if(msg) cat("\n####### Overall: ...  ", n.psi.ok, "selected breakpoint(s) \n\n")
    #cat(" ##### Overall ", nrow(r$psi),"selected breakpoints \n")
    #if(!is.list(r)) r<- olm
    bic.values=rbind(as.numeric(names(bicV)), bicV)
    rownames(bic.values)<-c("no. breakpoints", ICname)
    r$selection.psi <-list(bic.values=bic.values, npsi=n.psi.ok, cutvalues=cutvalues) ##cutvalues #including the extremes
    if(refit) r$selection.psi$psi.before.refit <-all.psi
    if(plot.ic) plot(t(bic.values), type="o"); points(n.psi.ok, min(bic.values[2,]), pch=19, cex=1.1)                 
    r
    #r<-segmented(olm, ~x, psi=r$psi[,"Est."], control=seg.control(n.boot=10,alpha=.01)) #fix.npsi=FALSE
    #r
  } #fine sel1()
  #=====================================================================
  # BIC<-function(obj, a=1){
  #   #Se a=1 questo e' il BIC classico (a meno di una costante)
  #   n <-length(obj$residuals)
  #   r <- n*log(sum(obj$residuals^2)/n) + (n-obj$df.residual)*(log(n)^a) #- 1
  #   r
  # }
  #=====================================================================
  type<-match.arg(type)
  if(!type%in%c("bic","aic") && Kmax!=2) stop("Kmax>2 is not (yet?) allowed with hypothesis testing procedures. Set type='bic' or 'aic'", call.=FALSE)
  if(!type%in%c("bic","aic") && G>1) stop("G>1 is allowed only with type='bic' or 'aic' ", call.=FALSE)
  #=====================================================================
  
  if(is.numeric(olm)){
    y<-olm
    if(missing(seg.Z)){
      Z<-x<- 1:length(y)
    } else {
      nomeX<-all.vars(seg.Z)
      if(length(nomeX)==1) {
        Z<- eval(parse(text=nomeX))
      } else {
        stop("a single segmented variable should be specified in 'seg.Z' ")
      }
    }
    olm <- lm(y~x)
  } else {
    if(!inherits(olm,"lm")) stop("'olm' does not appear a (g)lm fit")
    y<- model.response(model.frame(olm))
    if(is.matrix(y)){
      nomeY<-colnames(y)
    } else {
      nomeY<- all.vars(formula(olm))[1]
    }
    #browser()
    if(missing(seg.Z)){
      nomeX <- setdiff(all.vars(formula(olm)), nomeY)
      if(length(nomeX)==0 || length(nomeX)>1 || any(is.na(nomeX))) stop("I cannot determine the segmented variable")
      seg.Z<- as.formula(paste("~", nomeX ))
      Z <- olm$model[[nomeX]]
    } else {
      if(length(all.vars(seg.Z))>1) stop("Multiple variables are not allowed in seg.Z")
      nomeX<-all.vars(seg.Z)
      Z <- if(nomeX%in%all.vars(formula(olm))) olm$model[[nomeX]] else eval(parse(text=nomeX)) 
    }
  }
  m1 <-min(Z)
  m2 <-max(Z)

  if(is.null(th)) th <- diff(range(Z))/100

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
    if(type%in%c("bic","aic")){
      control1<-control
      control1$n.boot = 0
      control1$tol <- .001 #default e' .00001
      control1$alpha<-.01 #se lo aumenti puo' non funzionare bene se ci sono molti psi da selezionare..
      ICname<- if(type=="bic") "BIC" else "AIC"
      BIC.f<-if(type=="bic") BIC else AIC
      bicM0 <- BIC.f(olm)
      
      Kmax <- min(floor((olm$df.residual-1)/2), Kmax)
      
      npsi<-1:Kmax
      startpsi<-vector("list", length(npsi))  
      conv<-bic.values<- rep(NA, length(npsi))
      
      if(!is.null(olm$call$data)) assign(paste(olm$call$data), eval(olm$call$data, envir=parent.frame() ))
      
      npsiVeri<-0
      #fit with 1 breakpoint
      .a<-capture.output(os<- suppressWarnings(try(segmented(olm, seg.Z, npsi=1, control=control1), silent=TRUE)))
      ris<- NULL
      ris[[1]] <- os #<- suppressWarnings(try(segmented(olm, seg.Z, npsi=1, control=control1), silent=TRUE))
      #if fails try boot restating
      if(inherits(os, "try-error")) {
        .a <- capture.output(os<- suppressWarnings(try(segmented(olm, seg.Z, npsi=1, control=control), silent=TRUE)))
        ris[[1]]<- os
      }
      
      if(inherits(os, "segmented")){
        #if(is.null(th)) th <- drop(diff(os$rangeZ)/100)
        bic.values[1]<- BIC.f(os)
        Z<- os$model[,os$nameUV$Z]
        estpsi <- os$psi[,"Est."]
        M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
        psi0 <- sum(M[which.max(apply(M,1,diff)),])/2
        startpsi[[1]] <- sort(c(estpsi, psi0))
        conv[1] <- 1
        npsiVeri[length(npsiVeri)+1]<- length(estpsi)
      } else {
        estpsi <- (m1+m2)/2
        M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
        #psi0 <- sum(M[which.max(apply(M,1,diff)),])/2
        startpsi[[1]] <- estpsi
        bic.values[1]<- BIC.f(olm)+1
        conv[1] <- 0
        npsiVeri[length(npsiVeri)+1]<- 1
      }
      
      i=1 #ponilo =1 
      if(Kmax>=2){
        if(msg) {
          flush.console()
          cat(paste("No. of breakpoints: "))
        }
        earlyStop<- FALSE
        #==========================================================================================
        #= inizio for
        #browser()
        for(i in 2:Kmax){
          #source("C:/dati/lavori/segmented/segIntermedio/segmented/R/selgmented.R")
          #if(i==3) browser()
          .a <- capture.output(os<-suppressWarnings(try(segmented(olm, seg.Z, psi=startpsi[[i-1]], control=control1), silent=TRUE)))
          if(msg) {
            flush.console()
            cat(paste(i,".. "))
          }
          if(inherits(os, "segmented")) {
            conv[i]<-1
            estpsi <- os$psi[,"Est."]
            bic.values[i]<- BIC.f(os) #-2*logLik(ris[[i]]))+ edf*log(n)*Cn
            ris[[length(ris)+1]]<- os
            npsiVeri[length(npsiVeri)+1]<-length(estpsi)
            id<- which(diff(estpsi)<=th)+1
            if(length(id)>0){ #se ci sono psi troppo vicini..
              #elimina quelli "vicini" aggiungine altri in modo da stimare un altro modello con lo stesso numero di psi
              #Quindi alla fine dovresti avere 2 o piu' bic per uno stesso numero di breakpoints
              estpsi <- estpsi[-id]
              M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
              diffpsi <- apply(M,1,diff)
              psi0<-NULL
              for(j in 1:(length(id))) {
                psi0[length(psi0)+1] <- sum(M[which(diffpsi==rev(sort(diffpsi))[j]),])/2
              }
              startpsi[[i]] <- sort(c(estpsi, psi0))
            } else {
              #aggiungi uno starting psi
              M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
              diffpsi <- apply(M,1,diff)
              id.max.diffpsi <- which.max(diffpsi)
              psi0 <- sum(M[id.max.diffpsi,])/2
              startpsi[[i]] <- sort(c(estpsi, psi0))
            }
          } else {
            #Vuoi provare a ri-stimarlo con il boot rest?
            #se non e' arrivato a convergenza:
            #   1) prima prova a cambiare gli starting psi
            #   2) prova con il boot restart
            M<-matrix(c(m1 ,rep(startpsi[[i-1]], each=2), m2), 
                      ncol=2, nrow=length(startpsi[[i-1]])+1, byrow=TRUE)
            diffpsi <- apply(M,1,diff)
            psi0 <- sum(M[which.max(diffpsi),])/2
            start.ora<- sort(c(psi0, M[-c(which.min(diffpsi), nrow(M)),2]))
            startpsi[[i-1]] <- start.ora 
            .a <- capture.output(os<-suppressWarnings(try(segmented(olm, seg.Z, psi=start.ora, control=control1), silent=TRUE)))
            if(!inherits(os, "segmented")) { #vai con il boot restrat
              .a <- capture.output(os<-suppressWarnings(try(segmented(olm, seg.Z, psi=start.ora, control=control), silent=TRUE)))
            }
            if(inherits(os, "segmented")) {
              conv[i]<-1
              estpsi <- os$psi[,"Est."]
              bic.values[i]<- BIC.f(os) #-2*logLik(ris[[i]]))+ edf*log(n)*Cn
              ris[[length(ris)+1]]<- os
              npsiVeri[length(npsiVeri)+1]<-length(estpsi)
              #aggiungi uno starting psi
              M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
              diffpsi <- apply(M,1,diff)
              id.max.diffpsi <- which.max(diffpsi)
              psi0 <- sum(M[id.max.diffpsi,])/2
              startpsi[[i]] <- sort(c(estpsi, psi0))
            } else {
              #se dopo 3 tentativi non e' arrivato a conv comunque aggiungi uno starting psi,
              #  il bic e' bic.precedente+1 e ris[[i]] sara' NA
              #control1$alpha<-.005
              ris[[length(ris)+1]]<- NA
              conv[i]<-0
              bic.values[i]<- bic.values[i-1]+1
              M<-matrix(c(m1 ,rep(startpsi[[i-1]], each=2), m2), 
                      ncol=2, nrow=length(startpsi[[i-1]])+1, byrow=TRUE)
              diffpsi <- apply(M,1,diff)
              psi0 <- sum(M[which.max(diffpsi),])/2
            
              #aggiunge un psi ai precedenti
              startpsi[[i]] <- sort(c(startpsi[[i-1]], psi0))
              npsiVeri[length(npsiVeri)+1]<- length(startpsi[[i-1]])
            }
          }
          #if(i==8) browser()
          #un controllo su bic values.. fermarsi SE gli ultimi K sono NA oppure sono crescenti!!
          #bic.values<-bic.values[!is.na(bic.values)]
          #ind1 e' se gli ultimi modelli non sono arrivati a convergenza per cui il bic e' stato aumentato di 1 rispetto al precedente..
          #poiche' ci possono essere piu' bic per uno stesso numero di break, ne devi considerare solo uno (il minimo)
          #altrimenti -se i bic per uno stesso numero di break sono decrescenti- la valutazione del bic crescente e' sballata
          #if(i==5) browser()
          bicValuesTest<-tapply(na.omit(c(bicM0,bic.values)), npsiVeri, min)
          ind1<-(i>=stop.if && all(diff(rev(na.omit(bicValuesTest)))[1:3]==-1))
          ind2<-(i>=stop.if&& length(bicValuesTest)>stop.if && all(rev(na.omit(diff(bicValuesTest)))[1:stop.if]>0))
          if(ind1 || ind2) {earlyStop<-TRUE;break}
        } #end in for(2 in Kmax)
      
        #browser()
        #if(msg) cat(paste(" =>", length(psi.withCut)-length(bicVa), "unevaluated model(s) due to", stop.if, "increasing A/BIC value(s)..\n"))
        
        
        if(msg) cat("\n")  
      } else { #se Kmax=1
        earlyStop<-FALSE
      }
      #npsiVeri <-sapply(ris, function(.x) if(is.list(.x))nrow(.x$psi) else NA ) #[1:i]
      #npsiVeri <- npsiVeri[!is.na(npsiVeri)]

      #if(length(npsiVeri)!=max(npsiVeri)) {
      #  id.psi.repl<-which(diff(npsiVeri)==0)
      #  sapply(id.psi.repl, function(.x) which(npsiVeri==.x))
      #}
      
      bic.valuesOrig <- bic.values 
      #bic.values<- bic.values#[1:i] #bic.values[!is.na(bic.values)]
      bic.values <- c(bicM0, bic.values)
      bic.values <- bic.values[1:length(npsiVeri)]

      names(bic.values)<-npsiVeri
      n.psi.ok<- npsiVeri[which.min(bic.values)]

      #browser()
      
      if(n.psi.ok==0){
        m<-matrix(NA,1,1, dimnames=list(NULL, "Est."))
        olm$selection.psi<- list(bic.values=bic.values, npsi=n.psi.ok)
        olm$psi<-m
        if(msg){
          if(earlyStop) cat(paste("(search truncated at ", i, " breakpoints due to increasing values of ",  ICname ,") \n", sep=""))
          cat(paste("\n",ICname, " to detect no. of breakpoints:\n",sep=""))
          print(bic.values)
          cat(paste("\nNo. of selected breakpoints: ", n.psi.ok, " \n"))
        }
        if(return.fit) return(olm) else return(list(bic.values=bic.values, npsi=n.psi.ok))
      }
      
      #browser()
      
      if(Kmax==n.psi.ok && msg) warning(paste("The best",ICname, "value at the boundary. Increase 'Kmax'?"), call.=FALSE, immediate. = TRUE)
      id.best<- which.min(bic.values[-1])
      
      if(!return.fit) {
        r <- list(bic.values=bic.values, npsi=n.psi.ok)
        return(r) 
      }
      #browser()
      r<- r0 <- ris[[id.best]]
      
      #browser()
      
      rm.after.check <- 0
      all.psi<- r$psi[,"Est."]
      if(check.dslope){
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
          .a <- capture.output(r <- suppressWarnings(try(update(r0), silent=TRUE))) #senza boot
          if(!inherits(r,"segmented")){ #facciamo un altro tentativo..
            #r0$call$psi=all.psi
            #control$alpha <- .005
            r0$call$control<-control #con boot
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
        n.psi.ok<-length(all.psi) #in realta' gia' c'e' "n.psi.ok"
      }
      
      if(refit && length(all.psi)>0){
        r$call$psi<- all.psi
        #control$alpha <- .005
        r$call$control<-quote(control) #con boot
        r <- update(r)
      }

      if(plot.ic) {
        if(rm.after.check!=0){
          warning(paste("Some psi have been removed due to check on the slope difference; the", ICname, "plot could be misleading"), call.=FALSE)
        }
        plot(npsiVeri, bic.values, xlab=" No. of breakpoints", ylab=ICname, type="o")
        points(n.psi.ok, min(bic.values), pch=19, cex=1.1)
      }
      #browser()
      if(msg){
        if(earlyStop) cat(paste("(search truncated at ", i, " breakpoints due to ", stop.if, " increasing values of ",  ICname ,") \n", sep=""))
        cat(paste("\n",ICname, " to detect no. of breakpoints:\n",sep=""))
        print(bic.values)
        add.msg <- if(rm.after.check==0) " \n" else paste(" (", rm.after.check, " breakpoint(s) removed due to small slope diff)\n", sep="")
        cat(paste("\nNo. of selected breakpoints:", n.psi.ok, add.msg))
      }
      
      bic.values=rbind(npsiVeri, bic.values)
      rownames(bic.values)<-c("no. breakpoints", paste(ICname, "value",sep=""))
      r$selection.psi <- list(bic.values=bic.values, npsi=n.psi.ok)
      if(refit) r$selection.psi$psi.before.refit <-all.psi
      return(r)
      } else {
        #end aic/bic. Quindi se score o davies
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
                #olm$call$data<-quote(MF)
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
        n.psi.ok<-length(out$psi[,"Est."])
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
      }#end if(score o davies)
    } else { #se G>1
      r <- sel1(y=y, x=Z, G=G, Kmax=Kmax, type=type, th=th, refit=refit, check.dslope = check.dslope, 
              msg=msg, bonferroni=bonferroni, olm=olm, control=control)
    r
    }
}


