lines.segmented<-function(x, term, bottom=TRUE, shift=FALSE, conf.level=0.95, k=50, 
  pch=18, rev.sgn=FALSE,.vcov=NULL, .coef=NULL,...){
  if(missing(term)){
          if(length(x$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-x$nameUV$Z}
               }
  ss<-list(...)
  metodo<- if(!is.null(ss$method)) ss$method else "delta"
  colore<- if(is.null(ss$col)) 1 else ss$col
  usr <- par("usr")
  h<-(usr[4]-usr[3])/abs(k)
  y<- if(bottom) usr[3]+h else usr[4]-h
  m<- confint.segmented(object=x,parm=term,level=conf.level,rev.sgn=rev.sgn,digits=15,method=metodo,.vcov=NULL, .coef=NULL)
  #m<-r[[term]]
  #FORSE non e' necessaria
  #if(rev.sgn) m<- -m
  #ma invece serve il seguente (se length(psi)=1 e rev.sgn=T):
  m<-matrix(m,ncol=3)
  if(nrow(m)>1) m<-m[order(m[,1]),]
  est.psi<-m[,1]
  lower.psi<-m[,2]
  upper.psi<-m[,3]
  if(length(est.psi)>1) {
      y<- if(shift) y+seq(-h/2,h/2,length=length(est.psi)) else rep(y,length(est.psi))
      }
  #segments(lower.psi, y, upper.psi, y, ...)
  arrows(lower.psi, y, upper.psi, y, code=3, angle=90, length=.07, ...)
  points(est.psi,y,type="p",pch=pch,col=colore)
  }
