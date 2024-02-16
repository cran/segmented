pred.step.plot<-function(object, k=NULL, apprx=c("cdf","abs"), nomeZ, ...){
  apprx=match.arg(apprx)
  X=model.matrix.stepmented(object, k=k, apprx=apprx)
  V=vcov.stepmented(object)
  
  #browser()
  interc=TRUE
  nomiZ<- object$nameUV$Z
  nomiV<- object$nameUV$V
  nomiU<- object$nameUV$U
  nomiPsi<- gsub("V","psi", nomiV)
  #id.noV<-setdiff(colnames(X), nomiPsi)
  
  #se ci sono piu' variabili segmented seleziona solo i termini di un unica variabile
  if(missing(nomeZ)) nomeZ <-nomiZ[1]
  nomeU.ok   <- grep(paste(".", nomeZ,sep=""), nomiU, value=TRUE)
  nomePsi.ok <- grep(paste(".", nomeZ,sep=""), nomiPsi, value=TRUE)
  
  
  N=10
  vv<-matrix(, N, length(nomePsi.ok))
  for(i in 1:length(nomePsi.ok)){
    #psi <- object$psi[nomePsi.ok[i],"Est."]
    psi <- object$psi.rounded[1, nomePsi.ok[i]]
    se  <- object$psi[nomePsi.ok[i],"St.Err"]
    vv[,i] <-  qnorm(seq(.0005, .9995,l=N), psi, se)
  }
  vv<-as.vector(vv)  
  mi=object$rangeZ[1,nomeZ]
  ma=object$rangeZ[2, nomeZ]
  vv<-sort(c(object$psi.rounded[1, nomePsi.ok], seq(mi,ma,l=30), 
             vv))
  
  Xok<-matrix(, length(vv), 2*length(nomePsi.ok))
  colnames(Xok)<- c(nomeU.ok,nomePsi.ok)
  for(i in 1:length(nomePsi.ok)){
    Xok[, nomeU.ok[i]]   <- spline(object$Z[,nomeZ], X[, nomeU.ok[i]], xout=vv)$y
    Xok[, nomePsi.ok[i]] <- spline(object$Z[,nomeZ], X[, nomePsi.ok[i]], xout=vv)$y
  }
  
  id.ok = grep(paste(".", nomeZ,sep=""), colnames(X))
  
  psii = object$psi.rounded[1, nomePsi.ok]
  #psii = object$psi[nomePsi.ok,"Est."]
  M=sapply(psii, function(.x) 1*(vv>.x))
  cof<-object$coefficients[nomeU.ok]
  if(interc && "(Intercept)" %in% colnames(X)) {
    id.ok<-c(1, id.ok)
    Xok<-cbind(1, Xok)
    M = cbind(1, M)
    cof<-c(object$coefficients[1], cof)
  }
  
  se=sqrt(rowSums((Xok%*%V[id.ok,id.ok])*Xok))
  ff<-drop(M %*% cof)
  g <- cut(vv, breaks = c(mi, psii, ma), labels =FALSE, include.lowest = TRUE)
  r<-list(values=vv, fit=ff, se.fit=se, g=g)
  r
  #browser()
  #matplot(vv, cbind(ff, ff-2*se, ff+2*se), type="l", lty=c(1,2,2), col=1, xlab=colnames(object$Z)[1])
  #browser()
  #matpoints(vv, cbind(ff, ff-2*se, ff+2*se), type="l", lty=c(1,2,2), col=2)
  #matpoints(vv, cbind(ff, ff-2*se, ff+2*se), col=4, pch=4)
}

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
