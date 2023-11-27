points.segmented <-function(x, term, interc=TRUE, link=TRUE, 
    rev.sgn=FALSE, transf=I, .vcov=NULL, .coef=NULL, const=0, v=TRUE, ...){
#--------------
        f.U<-function(nomiU, term=NULL){
        #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
        #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
#-------------
      if(missing(term)){
          if(length(x$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-x$nameUV$Z}
               }
      opz<-list(...)
      if(is.null(opz$col)) opz$col <- 2
      if(is.null(opz$lty)) opz$lty <- 3
      nameV<- x$nameUV$V[f.U(x$nameUV$V, term)]
      psii<- x$psi[nameV, "Est."]
      d<-data.frame(a=psii)
      names(d)<-term
      opz$y<-broken.line(x,d, se.fit=FALSE, interc=interc, link=link, .coef=.coef, .vcov=.vcov)[[1]]
      #browser()
      if(rev.sgn) psii<- -psii
      opz$x<- psii 
      if(is.null(opz$cex)) opz$cex<-1.25
      if(is.null(opz$lwd)) opz$lwd<-1.6
      opz$y <- opz$y + const
      opz$y<-do.call(transf, list(opz$y))
      do.call(points, opz)
      if(v) segments(psii, par()$usr[3], psii, opz$y, lty = opz$lty, col=opz$col )
      invisible(NULL)
      }

      
      

      