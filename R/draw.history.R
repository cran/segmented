draw.history<-function(obj,term,...){
      if(missing(term)){
          if(length(obj$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-obj$nameUV$Z}
               }
      range.ok<-obj$rangeZ[,term]
      id.ok<-grep(term, rownames(obj$psi), extended=FALSE)
      est.psi<-obj$psi[id.ok,2]

      A<-matrix(unlist(obj$psi.history),nrow=nrow(obj$psi),byrow=FALSE)
      rownames(A)<-rownames(obj$psi)

      matplot(1:ncol(A),t(A)[,id.ok],type="b",pch=1:length(est.psi),col=1,
        xlab="iterations",ylab=term,ylim=range.ok,xaxt="n",...)
      axis(1,at=1:ncol(A))
      #if(rug) points(rep(1))
      abline(h=est.psi,lty=3)
      }
      

#da cancellare??!??!!





