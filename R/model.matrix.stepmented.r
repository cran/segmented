model.matrix.stepmented<-function(object, type=c("cdf","abs","none"), k=NULL, ...){ #ret.rangeZ=FALSE
  #if(!inherits(object, "segmented")) stop("A 'segmented' fit is requested")
  #browser()
  type=match.arg(type)
  if(type=="abs") stop(" type='abs' not (yet?) implemented")
  if(inherits(object, "lm")) {
    X<- qr.X(object$qr, ...)
    #if(inherits(object, "glm") ) {
    if(!is.null(object$weights)) { #questo vale sia per glm che per lm con weights  
      #W<-chol(diag(object$weights))
      #X <- X/diag(W)
      X<- X/sqrt(object$weights)
    }
  } else {
    class(object)<-class(object)[-1]
    X<-try(model.matrix(object,...), silent=TRUE)
    if(!is.matrix(X)) X<- model.matrix(object, data=model.frame(object))
  }
  p=ncol(X)
  n=nrow(X)
  nomiZ<- object$nameUV$Z
  nomiV<- object$nameUV$V
  nomiU<- object$nameUV$U
  nomiPsi<- gsub("V","psi", nomiV)
  id.noV<-setdiff(colnames(X), nomiPsi)
  
  #se object viene da stepmented() la matrice restituita e' sbagliata
  #se da stepreg() allora ok..
  dropV=FALSE
  if(strsplit(paste(object$call[[1]]), "\\.")[[1]][1]=="stepmented" && type=="none"){
    type<-"cdf"
    k=-500
    dropV=TRUE
  }
  
  if(type=="none") return(X[,id.noV, drop=FALSE]) #se object viene da stepmented() la matrice restituita e' sbagliata
  #se da stepreg() allora ok..
  
  if(inherits(object, "glm")) {
    sigma = if(object$family$family%in%c("poisson","binomial")) 1 else sqrt(object$deviance/object$df.residual)
  } else {
    sigma = sqrt(sum(object$residuals^2)/object$df.residual) 
  }
  
  #browser()
  
  maxZ.list<-NULL
  for(i in 1:length(nomiU)){
    nomeZ<- gsub("U[1-9]*[0-9].","",nomiU[i])
    Z<-object$Z[,nomeZ]
    minZ<-min(Z)
    maxZ<-max(Z)
    psi<-object$psi[nomiPsi[i],"Est."]

    if(type%in%c("cdf","abs" )){
      Z01<- (Z-minZ)/(maxZ-minZ)
      psi01<- (psi-minZ)/(maxZ-minZ)
      if(is.null(k)){
        idU<-match(nomiU[i],nomiU)
        snr.idU<-abs(object$coefficients[nomiU][idU])/sigma
        #ss01=n^(-(.6 + .3*log(snr.idU) -abs(psi01-.5)^(1/2)/n^(1/2)))
        ss01=n^(-(.6 + .5*log(snr.idU)/sqrt(snr.idU) -abs(psi01-.5)^(1/2)/n^(1/2)))
      } else {
        ss01=n^k
      }
      ss<- ss01*(maxZ-minZ)
      
      if(type=="cdf"){
        X[, nomiU[i]]<-  pnorm((Z-psi)/ss)
        if(nomiPsi[i]%in%colnames(X)) {
          X[, nomiPsi[i]] <- -(object$coefficients[nomiU[i]]/ss)*dnorm((Z-psi)/ss)
        } else {
          nomicolsX<-colnames(X)
          X <- cbind(X, -(object$coefficients[nomiU[i]]/ss)*dnorm((Z-psi)/ss))
          colnames(X)<- c(nomicolsX, nomiPsi[i] )
        }
        
      } else {
        xx <- Z-psi
        den <- -xx+2*xx*pnorm(xx/ss)+2*ss*dnorm(xx/ss) #.05*log(cosh((x-.5)/.05)))
        V <- (1/(2 * den))
        X[, nomiU[i]]<- (Z * V + 1/2)
        if(nomiPsi[i]%in%colnames(X)) {
          X[, nomiPsi[i]] <- -object$coefficients[nomiU[i]]*V
        } else {
          nomicolsX<-colnames(X)
          X <- cbind(X, -object$coefficients[nomiU[i]]*V)
          colnames(X)<- c(nomicolsX, nomiPsi[i] )
        }
      }
    }
    maxZ.list[[length(maxZ.list)+1]] <- maxZ-minZ
  }
  #browser()
  #if(ret.rangeZ) attr(X, "rangeZ")<- maxZ.list
  if(dropV) X<-X[, id.noV, drop=FALSE]
  return(X)
}

