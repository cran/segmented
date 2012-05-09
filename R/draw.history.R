draw.history<-function(obj,term,...){
#show.history() se c'è stato boot restart potrebbe produrre un grafico 2x1 di "dev vs it" and "no.of distinct vs it"
      if(missing(term)){
          if(length(obj$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-obj$nameUV$Z}
               }
      range.ok<-obj$rangeZ[,term]
      id.ok<-grep(paste("\\.",term,"$",sep=""),  rownames(obj$psi),value=FALSE)
      est.psi<-obj$psi[id.ok,2]

      A<-matrix(unlist(obj$psi.history),nrow=nrow(obj$psi),byrow=FALSE)
      rownames(A)<-rownames(obj$psi)

      matplot(1:ncol(A),t(A)[,id.ok],type="b",pch=1:length(est.psi),col=1,
        xlab="iterations",ylab=paste("breakpoint ","(",term,")",sep=""),
        ylim=range.ok,xaxt="n",...)
      axis(1,at=1:ncol(A),cex.axis=.7)
      #if(rug) points(rep(1))
      abline(h=est.psi,lty=3)
      }
      
