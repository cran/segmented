vcov.stepmented<-function(object, k=NULL, return.X=FALSE, zero.cor=TRUE, ...){
  #farlo su scala logit???
  #k=-1/2.. conservativi
  #k=-2/3.. ok, ma forse troppo piccoli..
  #================
  # U.ef<-U.pert<-matrix(, B, p)
  # if(!missing(seed)) set.seed(seed)
  # for(i in 1:B){
  #  id<-sample(n, replace=TRUE)
  #  U.ef[i,]<-colSums(U[id,]) #EF boot
  #  w<-rgamma(n,1,1) #o anche rnorm(n,1,1)
  #  U.pert[i,]<-colSums(U*w)
  # }
  # I1<-var(U.ef)
  # I2<-var(U.pert)
  # s2<-sum(o$residuals^2)/(n-p)
  # r<-list(teor= s2*crossprod(X), emp=crossprod(U), ef=I1,pert=I2)
  
  
  # var.Tay<-function(est1,est2,v1,v2,v12){
  #   r<- est1/est2
  #   vv<-(v1+v2*r^2-2*r*v12)/est2^2
  #   vv}
  # p<- length(object$coefficients)
  # n<- length(object$residuals)
  # objW<- object$objW
  # nomiCov<-names(object$coefficients)
  # nomiU <- object$nameUV$U
  # nomiVxb <- object$nameUV$V
  # nomiVxb<-gsub("V","psi",nomiVxb)
  # #browser()
  # 
  # if(!inherits(object, "lm")) stop(" A (g)lm object fit is requested")
  # if(inherits(object, "glm")){
  #   Cov<- chol2inv(object$obj.ok$qr$qr) #per GLM ma manca la dispersion..
  #   dispersion <- if(object$family$family%in%c("poisson","binomial")) 1 else object$dev/object$df.residual
  #   Cov<- Cov*dispersion
  # } else { #se lm
  #   w <- object$weights
  #   s2 <- if(is.null(w))  sum((object$residuals^2))/object$df.residual else sum(( w* object$residuals^2)[w > 0])/object$df.residual
  #   #Cov <- s2*summary.lm(objW)$cov.unscaled #chol2inv(qr.R(objW$qr))
  #   Cov<- s2 * solve(object$obj.ok$XtX)
  # }
  # Cov <- cbind(Cov,matrix(0, nrow=nrow(Cov), ncol=length(nomiVxb)))
  # Cov <- rbind(Cov,matrix(0, ncol=ncol(Cov), nrow=length(nomiVxb)))
  # rownames(Cov)<-colnames(Cov)<- nomiCov
  # diag(Cov)[nomiVxb]<- object$psi[,"St.Err"]^2
  # return(Cov)
  #================

  # fOb<-function(par, Xlin, Z, w=1){
  #   pLin<-ncol(Xlin)
  #   n<-nrow(Xlin)
  #   s<-1/sqrt(n)
  #   bLin <- par[1:pLin]
  #   muLin <-Xlin%*%bLin
  #   nPsi <- ncol(Z)
  #   diffSlope<-par[(pLin+1):(pLin+nPsi)]
  #   psi<-par[(pLin+nPsi+1):(pLin+2*nPsi)]
  #   Xpsi<-sapply(1:nPsi, function(.j){pnorm((Z[,.j] - psi[.j])/s)})
  #   muPsi<- Xpsi%*%diffSlope
  #   mu<-muLin+muPsi
  #   y<-object$fitted.values+object$residuals
  #   r<-sum(w*(y-mu)^2)/2
  #   r
  # }
  
 
  X0<-X<-model.matrix(object) #qr.X(object$qr) piu efficiente?
  p=ncol(X)
  n=nrow(X)
  nomiZ<- object$nameUV$Z
  nomiV<- object$nameUV$V
  nomiU<- object$nameUV$U
  nomiPsi<- gsub("V","psi", nomiV)
  id.noV<-setdiff(colnames(X), nomiPsi)
  
  #browser()
  #Xlin<-X[,setdiff(colnames(X),c(nomiU, nomiPsi))]
  #Hess<-optimHess(object$coefficients, fOb, Xlin=Xlin, Z=object$Z)
  ff<-function(k, se.only=TRUE, wght=FALSE, only.lin=FALSE, return.X=FALSE){
    #browser()
    maxZ.list<-NULL
    for(i in 1:length(nomiU)){
      nomeZ<- gsub("U[1-9].","",nomiU[i])
      Z<-object$Z[,nomeZ]
      minZ<-min(Z)
      maxZ<-max(Z)
      psi<-object$psi[nomiPsi[i],"Est."]
      Z<- (Z-minZ)/(maxZ-minZ)
      psi<- (psi-minZ)/(maxZ-minZ)
      ss=n^k
      if(wght){
        idU<-match(nomiU[i],nomiU)
        snr.idU<-abs(object$coefficients[nomiU][idU])/sigma
        ss=n^(-(.6 + .3* log(snr.idU) -abs(psi-.5)^.5/sqrt(n)))
        #ss=n^(-(.6 + .07* log(snr.idU)*log10(n) -abs(psi-.5)^.5/sqrt(n)))
        #.6 + .3* log(o$coefficients[2]/s) -abs(o$psi[,1]-.5)^.5/sqrt(n) +log(log(log10(n)))/3
      }
      #browser()
      X0[,nomiU[i]]<- 1*(Z>psi)
      X[,nomiU[i]]<- pnorm((Z-psi)/ss)
      X0[, nomiPsi[i]] <- X[, nomiPsi[i]] <- -(object$coefficients[nomiU[i]]/ss)*dnorm((Z-psi)/ss)
      maxZ.list[[length(maxZ.list)+1]]<-maxZ-minZ
    }
    if(return.X) return(list(X0,X))
    #X<-X[, -match(nomiPsi, colnames(X))]
    #grad<-crossprod(X, obj)
    if(!inherits(object, "lm")) stop(" A (g)lm object fit is requested")
    if(inherits(object, "glm")){
      variance = object$family$variance
      linkinv = object$family$linkinv
      mu.eta = object$family$mu.eta
      
      eta <- object$linear.predictors
      prior.weights=object$prior.weights
      
      mu = linkinv(eta) 
      varg = variance(mu)
      invgprime = mu.eta(eta)
      ww<- prior.weights*(invgprime^2/varg) #object$prior.weights*object$weights
      if(only.lin){
        invH<- solve(crossprod(sqrt(ww)*X0[,id.noV,drop=FALSE]))
        U<-X0[,id.noV,drop=FALSE]*(prior.weights*object$residuals*invgprime/varg)
      } else {
        invH<- solve(crossprod(sqrt(ww)*X))
        U<-X0*(prior.weights*object$residuals*invgprime/varg)
      }
    } else { #se lm
      if(is.null(object$weights)){
        w<-1
      } else {
        w<-object$weights
      }
      if(only.lin){
        invH<- solve(crossprod(sqrt(w)*X0[,id.noV,drop=FALSE]))
        U<-X0[,id.noV,drop=FALSE]*(w*object$residuals)
      } else {
        invH<- solve(crossprod(sqrt(w)*X))
        U<-X0*(w*object$residuals)
      }
    }
    INF<- crossprod(U)
    V =invH %*% INF %*% invH
    if(only.lin) return(V)
    #browser()
    for(i in 1:length(nomiPsi)){
      V[,nomiPsi[i]]<-V[,nomiPsi[i]]*maxZ.list[[i]]
      V[nomiPsi[i],]<-V[nomiPsi[i],]*maxZ.list[[i]]
    }
    if(se.only) {
      V<-diag(V)[nomiU]
      V<-sum((V-seb)^2)
    }
  return(V)
  }
  #s2<-sum(object$residuals^2)/(n-p) #sqrt(diag(s2*invH))  funziona benino..
  #  invH1<- solve(Hess)
  #V1=invH1 %*% INF %*% invH1
  #V=list(sqrt(diag(V)), sqrt(diag(V1)))
  #V<-list(sqrt(diag(V)), sqrt(diag(s2*invH)))
 #se<-sqrt(diag(X%*%V%*%t(X)))
  #M<-cbind(ff-2*se, ff, ff+2*se)
  #browser()
  if(is.numeric(k)){
    if(length(k)==1){
      #if(!is.null(object$vcov) && k==-2/3) return(object$vcov)
      V<-ff(k, se.only=FALSE,return.X=return.X)
      } else {
        seb<- summary(object)$sigma*sqrt(diag(solve(object$obj.ok$XtX))[nomiU])
        .o<-optimize(ff, c(k[1], k[2]))
        V<-ff(.o$minimum, se.only=FALSE, return.X=return.X)
        attr(V, "k")<-.o$minimum
      }
  } else {
    #if(!is.null(object$vcov)) return(object$vcov)
    #seb<- summary(object)$sigma*sqrt(diag(solve(object$obj.ok$XtX))[nomiU])
    seb=0 #non serve in realta'...
    sigma= sqrt(sum(object$residuals^2)/object$df.residual)
    V<-ff(1, se.only=FALSE, wght=TRUE, return.X=return.X)
  }
  if(return.X) return(V)
  V0<-ff(1, only.lin=TRUE)
  #browser()
  V[id.noV,id.noV]<-V0
  if(zero.cor) V[nomiPsi, id.noV]<- V[id.noV, nomiPsi] <-0
  V 
}


  


