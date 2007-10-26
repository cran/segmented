plot.segmented<-function(x,term=NULL,se=FALSE,const=coef(x)["(Intercept)"], add=FALSE,
  linkinv=FALSE, h=300, show.gap=TRUE, ...){
#Questa è quasi uguale a ff nel file bp.r
#Non ha l'argomento y=NULL della generica plot (mhmm..anche plot.density() and plot.gam() non
#hanno l'argomento y  bhuu???)
    if(se) {
      se<-FALSE
      warning("se=TRUE not (yet) implemented",call. = FALSE)
      }
    mio.plot<-function(...,col,lty,lwd) plot.default(...)
    mio.lines<-function(...,ylim,xlim,pch,xlab,ylab) lines(...)
      if(missing(term)){
          if(length(x$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-x$nameUV$Z}
               }

      id.Z<-match(x$nameUV$Z,names(coef(x)))
      id.U<-match(x$nameUV$U,names(coef(x)))
      id.V<-match(x$nameUV$V,names(coef(x)))

      id.term<-grep(term, names(coef(x)), extended=FALSE) #tutti i coeff di "term"
      id.U<-id.term[id.term%in%id.U]
      id.V<-id.term[id.term%in%id.V]
      id.Z<-id.term[id.term%in%id.Z]

      cof<-c(coef(x)[id.U], coef(x)[id.V]*coef(x)[id.U])
      cof<- if(is.na(id.Z)||length(id.Z)==0) c(0,cof) else c(coef(x)[id.Z],cof)

      if(se){
      #oppure: (per aggiungere gli SE)
      A<-diag(c(1,1,coef(x)[id.U]))
      id.ok<-c(id.Z,id.U,id.V)
      cof<-A%*%coef(x)[id.ok]
      var.cof<-A%*%vcov(x)[id.ok,id.ok]%*%A
      }

      idpsi<-grep(term, rownames(x$psi), extended=FALSE)
      psi<-x$psi[idpsi,2]

      xx<-seq(min(x$rangeZ[,term]), max(x$rangeZ[,term]), length=h)
      nn<-length(xx)
      k<-length(psi)
      Z<-matrix(rep(xx,times=k),nrow=nn,byrow = FALSE)
      PSI<- matrix(rep(psi, rep(nn,k)), ncol=k)
      U<-pmax((Z -PSI), 0)
      V<-ifelse(Z>PSI,-1,0)
      X<-cbind(xx,U,V)
      yhat<-drop(X%*%cof)
      if(!show.gap){
        X <- cbind(xx, U)
        cof<-lm.fit(x=X,y=yhat)$coefficients
        yhat <- drop(X %*% cof)
        }
      if(se){
      #per la seguente espressione c'è una formula più efficiente..
        se.yhat<-sqrt(diag(X%*%tcrossprod(var.cof,X))) #X%*%var.cof%*%t(X)
        inf.yhat<-yhat-1.96*se.yhat
        sup.yhat<-yhat+1.96*se.yhat
        } else {
            yhat<-yhat+const}
      ylab<-"link(Fitted Values)"
      if(inherits(x,what="glm",which=FALSE) && linkinv) {
            yhat<-x$family$linkinv(yhat)
            if(se){
                inf.yhat<-x$family$linkinv(inf.yhat)
                sup.yhat<-x$family$linkinv(sup.yhat)
                }
            ylab<-"Fitted Values"
            }

      sel.g<-rowSums(V) #questa variabili definisce i gruppi (segmenti) di appartenenza
      #if(!add) plot.default(xx, yhat, type="n", ylab=ylab, xlab=term, ...)
      #for(i in sel.g) lines(xx[sel.g==i], yhat[sel.g==i])
      if(!add) {
        if(!se) {mio.plot(xx, yhat, type="n", ylab=ylab, xlab=term,...)}
            else {matplot(xx,cbind(yhat,inf.yhat,sup.yhat), type="n", ylab=ylab, xlab=term,...)}
        }
      for(i in sel.g) mio.lines(xx[sel.g==i], yhat[sel.g==i],...)
      if(se) {
        for(i in sel.g) mio.lines(xx[sel.g==i], inf.yhat[sel.g==i],...)
        for(i in sel.g) mio.lines(xx[sel.g==i], sup.yhat[sel.g==i],...)
        }
      
      }#end plot.segmented
