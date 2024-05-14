vcov.stepmented<-function(object, k=NULL, zero.cor=TRUE, type=c("cdf", "none", "abs"), ...){
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
  
  #browser()
  
  if(!inherits(object, "lm")) stop("A stepmented (g)lm is requested")
  #calcola il parametro di disp che serve sempre..
  if(inherits(object, "glm")){
    disp <- object$deviance/object$df.residual
    if(object$family$family %in% c("binomial","poisson")) disp <-1 
  } else {
    ww <- if(is.null(object$weights)) 1 else object$weights 
    disp <- sum(ww*object$residuals^2)/object$df.residual
  }
  #browser()
  
  type<-match.arg(type)
  if(type=="abs") stop("type='abs' not (yet?) implemented")
  b<-coef.stepmented(object, FALSE)
  b <- b[b!=0]
  pLin<-length(b)
  
  if(type=="none"){
    
    V<-if(is.null(object$obj.ok)) chol2inv(object$qr$qr[1:pLin, 1:pLin, drop = FALSE]) else object$obj.ok$invXtX
    #ok anche in presenza di pesi va bene!
    V<- V*disp
    #browser()
    colnames(V)<-rownames(V)<- names(b)
    return(V)
  }
  #=====================================================
  #X0 deve avere le variabili I(x>psi). Se object e' restituito da segreg allora model.matrix.stepmented(object, apprx = "no")
  #funziona, se restituito da stepmented no. Allora e' meglio mettere
  #X0<- model.matrix.stepmented(object, apprx = "no") #funziona solo con oggetti stepreg
  #browser()
  
  X0<- model.matrix.stepmented(object, type = "cdf", k=-100)[,1:pLin]
  X <- model.matrix.stepmented(object, type = type, k=k) #qr.X(object$qr) piu efficiente?
  
  X0 <- cbind(X0, X[,setdiff(1:ncol(X),1:ncol(X0))]) #aggiungi i termini relativi ai psi
  maxZ.list <- attr(X, "rangeZ") 
  attr(X, "rangeZ")<-NULL
  
  #browser()
  
  p=ncol(X)
  n=nrow(X)
  nomiZ<- object$nameUV$Z
  nomiV<- object$nameUV$V
  nomiU<- object$nameUV$U
  nomiPsi<- sub("V","psi", nomiV)
  id.noV<-setdiff(colnames(X), nomiPsi)
  
  #====================================
  #inutile trasf.X(), il lavoro lo fa model.matrix.stepmented
  # trasf.X<-function(k){
  #   #trasforma la matrice del disegno di un modello stepmented.. utile per il calcolo della vcov
  #   maxZ.list<-NULL
  #   for(i in 1:length(nomiU)){
  #     nomeZ<- gsub("U[1-9].","",nomiU[i])
  #     Z<-object$Z[,nomeZ]
  #     minZ<-min(Z)
  #     maxZ<-max(Z)
  #     psi<-object$psi[nomiPsi[i],"Est."]
  #     Z<- (Z-minZ)/(maxZ-minZ)
  #     psi<- (psi-minZ)/(maxZ-minZ)
  #     if(is.null(k)){
  #       idU<-match(nomiU[i],nomiU)
  #       snr.idU<-abs(object$coefficients[nomiU][idU])/sigma
  #       ss=n^(-(.6 + .3* log(snr.idU) -abs(psi-.5)^.5/sqrt(n)))
  #       #ss=n^(-(.6 + .07* log(snr.idU)*log10(n) -abs(psi-.5)^.5/sqrt(n)))
  #       #.6 + .3* log(o$coefficients[2]/s) -abs(o$psi[,1]-.5)^.5/sqrt(n) +log(log(log10(n)))/3
  #     } else {
  #       ss=n^k
  #     }
  #     #browser()
  #     X0[,nomiU[i]]<- 1*(Z>psi)
  #     X[,nomiU[i]]<- pnorm((Z-psi)/ss)
  #     X0[, nomiPsi[i]] <- X[, nomiPsi[i]] <- -(object$coefficients[nomiU[i]]/ss)*dnorm((Z-psi)/ss)
  #     maxZ.list[[length(maxZ.list)+1]]<-maxZ-minZ
  #   }
  #   return(list(X0=X0, X=X, maxZ.list=maxZ.list))
  # }
  #====================================
  #R<-trasf.X(k)
  #browser()
  
  #X0<-R$X0
  #X<- R$X
  
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
      # if(only.lin){
      #   invH<- solve(crossprod(sqrt(ww)*X0[,id.noV,drop=FALSE]))
      #   U<-X0[,id.noV,drop=FALSE]*(prior.weights*object$residuals*invgprime/varg)
      # } 
      invH<- solve(crossprod(sqrt(ww)*X))
      U<-X0*(prior.weights*object$residuals*invgprime/varg)
    } else {#se lm
      w <- if(is.null(object$weights) || sd(object$weights)==0) 1 else object$weights 
      # if(only.lin){
      #   invH<- solve(crossprod(sqrt(w)*X0[,id.noV,drop=FALSE]))
      #   U<-X0[,id.noV,drop=FALSE]*(w*object$residuals)
      # } 
      #browser()
      invH<- solve(crossprod(sqrt(w)*X))
      U<-X0*(w*object$residuals)
    }
  INF<- crossprod(U)
  V =invH %*% INF %*% invH
  #browser()
  # for(i in 1:length(nomiPsi)){
  #     V[,nomiPsi[i]]<-V[,nomiPsi[i]]*maxZ.list[[i]]
  #     V[nomiPsi[i],]<-V[nomiPsi[i],]*maxZ.list[[i]]
  # }

  if(zero.cor) V[nomiPsi, id.noV]<- V[id.noV, nomiPsi] <-0
  V 
}


  


