selgmented<-function(olm, seg.Z, alpha=0.05, type=c("score" , "davies", "bic", "aic"), control=seg.control(), refit=TRUE, 
            stop.if=6, return.fit=TRUE, bonferroni=FALSE, Kmax=2, msg=TRUE, plot.ic=FALSE, th=NULL){
  #refit.last: if TRUE (and return.fit=TRUE), the last (selected) model is re-fitted using the 'control' controlling options (e.g. with n.boot>0)
  #stop.if. If type="bic" or "aic", the search of number of breakcan be stopped when the last 'stop.if' fits provide higher aic/bic value
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
  
  if(is.numeric(olm)){
    y<-olm
    x<-1:length(y)
    olm<-lm(y~x)
  }
  if(missing(seg.Z)){
    nomeX <- all.vars(formula(olm))[2]
    if(length(nomeX)>1 || any(is.na(nomeX))) stop("I cannot determine the segmented variable")
    seg.Z<- as.formula(paste("~", nomeX ))
  } else {
    if(length(all.vars(seg.Z))>1) stop("Multiple variables are not allowed in seg.Z")
  }
  type<-match.arg(type)
  if(!type%in%c("bic","aic") && Kmax!=2) stop("Kmax>2 is not (yet?) allowed with hypothesis testing procedures", call.=FALSE)

  if(type%in%c("bic","aic")){
    ICname<- if(type=="bic") "BIC" else "AIC"
    control1<-control
    control1$n.boot=0
    control1$tol <- .001 #default e' .00001
    control1$alpha<-.01
    BIC.f<-if(type=="bic") BIC else AIC
    bicM0 <- BIC(olm)
    npsi<-1:Kmax
    startpsi<-ris<-vector("list", length(npsi))  
    conv<-bic.values<- rep(NA, length(npsi))
    
    #fit with 1 breakpoint
    ris[[1]]<- os<- suppressWarnings(try(segmented(olm, seg.Z, npsi=1, control=control1, silent=TRUE)))
    
    if(is.null(th)) th <- drop(diff(os$rangeZ)/100)
    
    #browser()

    bic.values[1]<- BIC.f(os)
    Z<- os$model[,os$nameUV$Z]
    m1 <-min(Z)
    m2 <-max(Z)
    estpsi <- os$psi[,"Est."]
    M<-matrix(c(m1 ,rep(estpsi, each=2), m2), ncol=2, nrow=length(estpsi)+1, byrow=TRUE)
    psi0 <- sum(M[which.max(apply(M,1,diff)),])/2
    startpsi[[1]] <- sort(c(estpsi, psi0))
    if(msg) {
      flush.console()
      cat(paste("No. of breakpoints: "))
    }
    #browser()
    for(i in 2:Kmax){
      #if(i==7) browser()
      ris[[i]]<-os<-suppressWarnings(try(segmented(olm, seg.Z, psi=startpsi[[i-1]], control=control1), silent=TRUE))
      if(inherits(os, "segmented")) {
        conv[i]<-1
        if(msg) {
          flush.console()
          cat(paste(i,".. "))
        }
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
    }
    #browser()
    #npsiVeri <- unlist(sapply(ris, function(x)nrow(x$psi))) #n.breakpoints nei diversi modelli
    npsiVeri <-c(1,unlist(sapply(startpsi, length))[1:(i-1)])
    
    bic.valuesOrig <- bic.values 
    bic.values<- bic.values[1:i] #bic.values[!is.na(bic.values)]
    bic.values <- c(bicM0, bic.values)
    names(bic.values)<-c("0", npsiVeri)
    n.psi.ok<- c(0,npsiVeri)[which.min(bic.values)]
    
    if(msg){
      cat(paste("\n\nNo. of selected breakpoints: ", n.psi.ok, " \n"))
      if(any(is.na(bic.valuesOrig))) cat(paste("(search truncated at ", max(npsiVeri), " breakpoints due to increasing values of ",  ICname ,") \n", sep=""))
      cat(paste("\n",ICname, " to detect no. of breakpoints:\n",sep=""))
      print(bic.values)
    }
    
    if(n.psi.ok==0){
      m<-matrix(NA,1,1, dimnames=list(NULL, "Est."))
      olm$psi<-m
      return(olm)
    }

    if(i==n.psi.ok) warning(paste("The best",ICname, "value at the boundary. Increase 'Kmax'?"), call.=FALSE)
    
    r<-list(bic.values=bic.values, n.psi=n.psi.ok)

    id.best<- which.min(bic.values[-1])
    
    #browser()
    
    if(!return.fit) {
      return(r) 
    } else {
    #ATTENZIONE: SE NON HA SELEZIONATO BREAKPOINTS???
      if(refit){
        ris[[id.best]]$call$psi<- ris[[id.best]]$psi[,"Est."]#startpsi[[n.psi.ok-1]]
        control$alpha <- .01
        ris[[id.best]]$call$control<-control
        ris[[id.best]] <- update(ris[[id.best]])
      }
      ris[[id.best]]$selection.psi <- bic.values #era ris[[n.psi.ok]]$selection.psi <- bic.values. cambiato nella 1.6-3
      if(plot.ic) {
        #bic.values<-bic.values[!is.na(bic.values)]
        #plot(0: (length(bic.values)-1), bic.values, xlab=" No. of breakpoints", ylab=ICname, type="o")
        plot(c(0,npsiVeri), bic.values, xlab=" No. of breakpoints", ylab=ICname, type="o")
        points(n.psi.ok, min(bic.values), pch=19)
      }
    return(ris[[id.best]])
    }
  }
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
}
