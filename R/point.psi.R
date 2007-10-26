point.psi<-function(obj,term,bottom=TRUE,conf.level=0.95,k=50,pch=18,rev.sgn=FALSE,...){
  if(missing(term)){
          if(length(obj$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-obj$nameUV$Z}
               }
  ss<-list(...)
  colore<- if(is.null(ss$col)) 1 else ss$col
  usr <- par("usr")
  h<-(usr[4]-usr[3])/abs(k)
  y<- if(bottom) usr[3]+h else usr[4]-h
  r<- ci.psi(obj,conf.level=conf.level,digits=15)
  m<-r[[term]]
  if(rev.sgn) m<- -m
  est.psi<-m[,1]
  lower.psi<-m[,2]
  upper.psi<-m[,3]
  segments(lower.psi, y, upper.psi, y, ...)
  points(est.psi,rep(y,length(est.psi)),type="p",pch=pch,col=colore)
  }
  
