model.matrix.stepmented<-function(object, k=NULL, apprx=c("cdf","abs"), ...){
  #if(!inherits(object, "segmented")) stop("A 'segmented' fit is requested")
  apprx=match.arg(apprx) 
  if(inherits(object, "lm")) {
    X<- qr.X(object$qr, ...)
    if(inherits(object, "glm")) {
      W<-chol(diag(object$weights))
      X <- X/diag(W)
    }
  } else {
    class(object)<-class(object)[-1]
    X<-try(model.matrix(object,...), silent=TRUE)
    if(!is.matrix(X)) X<- model.matrix(object, data=model.frame(object))
  }
  X0<-X
  p=ncol(X)
  n=nrow(X)
  nomiZ<- object$nameUV$Z
  nomiV<- object$nameUV$V
  nomiU<- object$nameUV$U
  nomiPsi<- gsub("V","psi", nomiV)
  id.noV<-setdiff(colnames(X), nomiPsi)
  maxZ.list<-NULL
  #browser()
  sigma=sqrt(sum(object$residuals^2)/object$df.residual)
  
  for(i in 1:length(nomiU)){
    nomeZ<- gsub("U[1-9].","",nomiU[i])
    Z<-object$Z[,nomeZ]
    minZ<-min(Z)
    maxZ<-max(Z)
    psi<-object$psi[nomiPsi[i],"Est."]
    Z01<- (Z-minZ)/(maxZ-minZ)
    psi01<- (psi-minZ)/(maxZ-minZ)
    if(is.null(k)){
      idU<-match(nomiU[i],nomiU)
      snr.idU<-abs(object$coefficients[nomiU][idU])/sigma
      ss01=n^(-(.6 + .3* log(snr.idU) -abs(psi01-.5)^.5/sqrt(n)))
    } else {
      ss01=n^k
    }
    ss<- ss01*(maxZ-minZ)
    
    X0[, nomiU[i]]<-  pnorm((Z-psi)/ss) 
    X0[, nomiPsi[i]] <- -(object$coefficients[nomiU[i]]/ss)*dnorm((Z-psi)/ss)
    
    #X0[, nomiU[i]]<-  pnorm((Z01-psi01)/ss01) #1*(Z>psi)
    #X0[, nomiPsi[i]] <- -(object$coefficients[nomiU[i]]/ss01)*dnorm((Z01-psi01)/ss01)
    
    #opzione 2:
    xx <- Z-psi
    den <- -xx+2*xx*pnorm(xx/ss)+2*ss*dnorm(xx/ss) #.05*log(cosh((x-.5)/.05)))
    #den <- abs(xx)
    #browser()
    V <- (1/(2 * den))
    
    X[,nomiU[i]]<- (Z * V + 1/2) #U <-
    X[, nomiPsi[i]] <- -object$coefficients[nomiU[i]]*V
    maxZ.list[[length(maxZ.list)+1]] <- maxZ-minZ
  }
  #browser()
  #X e' basata sull'approx del valore assoluto
  XX<-if(apprx=="cdf") X0 else X
  return(XX)
}

