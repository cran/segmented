#new predict.segmented
predict.segmented<-function(object, newdata, se.fit=FALSE, interval=c("none","confidence", "prediction"), 
            type = c("link", "response"),# "terms"),
            level=0.95, .coef=NULL, ...){
  blockdiag <- function(...) {
    args <- list(...)
    nc <- sapply(args,ncol)
    cumnc <- cumsum(nc)
    ##  nr <- sapply(args,nrow)
    ## NR <- sum(nr)
    NC <- sum(nc)
    rowfun <- function(m,zbefore,zafter) {
      cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
            matrix(0,ncol=zafter,nrow=nrow(m)))
    }
    ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
    for (i in 2:length(args)) {
      ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
    }
    ret
  }
  
  dummy.matrix<-function(x.values, x.name, obj.seg, psi.est=TRUE, isV=FALSE, .coef=NULL){ 
    #given the segmented fit 'obj.seg' and a segmented variable x.name with corresponding values x.values,
    #this function simply returns a matrix with columns (x, (x-psi)_+, -b*I(x>psi))
    #or  ((x-psi)_+, -b*I(x>psi)) if obj.seg does not include the coef for the linear "x"
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
    estcoef <- if(is.null(.coef)) coef(obj.seg) else .coef
    if(length(isV)==1) isV<-c(FALSE,isV)
    n<-length(x.values)
    #le seguenti righe selezionavano (ERRONEAMENTE) sia "U1.x" sia "U1.neg.x" (se "x" e "neg.x" erano segmented covariates)
    #nameU<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$U, value = TRUE)
    #nameV<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$V, value = TRUE)
    nameU<-obj.seg$nameUV$U[f.U(obj.seg$nameUV$U,x.name)]
    nameV<-obj.seg$nameUV$V[f.U(obj.seg$nameUV$V,x.name)] #grep(x.name, obj.seg$nameUV$V, value = TRUE)
    
    #browser()
    if(is.null(obj.seg$constr)){
      diffSlope<-estcoef[nameU]
    } else {
      diffSlope<-drop(obj.seg$constr$invA.RList[[match(x.name, obj.seg$nameUV$Z)]]%*%estcoef[nameU])[-1]
    }
    
    est.psi<-obj.seg$psi[nameV,"Est."]
    se.psi<-obj.seg$psi[nameV, "St.Err"]
    k<-length(est.psi)
    PSI <- matrix(rep(est.psi, rep(n, k)), ncol = k)
    SE.PSI <- matrix(rep(se.psi, rep(n, k)), ncol = k)
    newZ<-matrix(x.values, nrow=n,ncol=k, byrow = FALSE)
    
    
    dummy1<-if(isV[1]) (newZ-PSI)*pnorm((newZ-PSI)/SE.PSI) else  (newZ-PSI)*(newZ>PSI) #pmax(newZ-PSI,0)
    
    if(psi.est){
      V<-if(isV[2]) -pnorm((newZ-PSI)/SE.PSI) else -(newZ>PSI) #ifelse(newZ>PSI,-1,0)
      dummy2<- if(k==1) V*diffSlope  else V%*%diag(diffSlope) #t(diffSlope*t(-I(newZ>PSI)))
      newd<-cbind(x.values,dummy1,dummy2)
      #colnames(newd)[1]<- x.name
      colnames(newd)<-c(x.name,sub("psi","U", nameV), nameV)
    } else {
      newd<-cbind(x.values,dummy1)
      #colnames(newd)[1]<- x.name
      colnames(newd)<-c(x.name, sub("psi","U", nameV))
    }
    
    #if(!x.name%in%names(coef(obj.seg))) newd<-newd[,-1,drop=FALSE] #restituisce sempre il termine principale..
    #aggiungi (eventualmente) le colonne relative ai psi noti
    all.psi<-obj.seg$indexU[[x.name]]
    if(length(all.psi)!=k){
      nomi.psi.noti<-setdiff(names(all.psi),nameU)
      psi.noti<-setdiff(all.psi, est.psi)
      PSI.noti <- matrix(rep(psi.noti, rep(n, length(psi.noti))), ncol = length(psi.noti))
      nomi<-c(colnames(newd),nomi.psi.noti)
      newd<-cbind(newd, (newZ-PSI.noti)*(newZ>PSI.noti))
      colnames(newd)<-nomi
    }
    return(newd)
  }
  estcoef <- if(is.null(.coef)) coef(object) else .coef
  nLin<- length(setdiff(names(coef(object)), c(object$nameUV$U,object$nameUV$V)))
  nSeg<- length(object$nameUV$Z)
  type<-match.arg(type)
  interval<-match.arg(interval)
  #browser()
  if(inherits(object, "glm") && object$family$family!="gaussian" && interval=="prediction") 
    stop("prediction intervals are not allowed with non-gaussian glm")
  nameU<-object$nameUV$U
  nameV<-object$nameUV$V
  nameZ<-object$nameUV$Z
  
  
  
  if(missing(newdata)){
    X <- model.matrix(object)
  } else {
    #browser()
    nomiLin <- setdiff(all.vars(formula(object))[-1], c(object$nameUV$U,object$nameUV$V))
    if(any(is.na(match(nomiLin, names(newdata))))) stop(" 'newdata' should includes all variables")
    #devi trasformare la variabili segmented attraverso dummy.matrix()
    #browser()
    n<-nrow(newdata)
    r<-NULL
    for(i in 1:length(nameZ)){
      x.values<-newdata[[nameZ[i]]]
      DM<-dummy.matrix(x.values, nameZ[i], object)
      #dummy.matrix() non restituisce i nomi.. quindi devi aggiungerli
      #nomiV.i<- grep(nameZ[i], object$nameUV$V, value = TRUE)
      #nomiU.i<- sub("psi","U", nomiV.i)
      #nomiOK <- c(nomiU.i,nomiV.i)
      #colnames(DM)<- if(nameZ[i]%in%colnames(DM)) c(nameZ[i])
      r[[i]]<-DM
    }
    newd.ok<-data.frame(matrix(unlist(r), nrow=n, byrow = FALSE))
    names(newd.ok)<- unlist(sapply(r, colnames))
    idZ<-match(nameZ, names( newdata))
    X<-data.matrix(cbind(newdata[,-idZ, drop=FALSE], newd.ok)) 
    if("(Intercept)" %in% names(estcoef)) X<-cbind("(Intercept)"=1,X)
  }
  
  #browser()
  colnomi<- colnames(X)

  if(!is.null(object$constr)){
    for(i in 1:length(nameZ)){
        nomeU.i<-grep(object$nameUV$Z[i], object$nameUV$U, value=TRUE)
        idU.i <- match(nomeU.i, names(estcoef))
        coef.new<-drop(object$constr$invA.RList[[i]]%*%estcoef[nomeU.i])
        names(coef.new)<-c(object$nameUV$Z[i], 
                           paste("U",1:(length(coef.new)-1),".",object$nameUV$Z[i],sep="" ))
        estcoef<-append(estcoef[-idU.i], coef.new, after=idU.i[1]-1)
    }
  } else {
    for(i in 1:length(nameZ)){
      if(!nameZ[i]%in%names(estcoef)) colnomi<-setdiff(colnomi, nameZ[i])
    }
  }
  
  colnomi.noV <- setdiff(colnomi, nameV)
  X.noV <- X[, colnomi.noV]
  
  #browser()
  estcoef.noV<- estcoef[colnomi.noV]
  #estcoef.noV<- estcoef[-match(nameV,names(estcoef), 0)]
  mu <- eta<- drop(X.noV%*% estcoef.noV)
  
  #ATTENZIONE c'e' il problema dell'appaiamento dei nomi!!!
  #il problema e' che estcoef non ha sempre nomi!! 

  X <- X[,c(colnomi.noV, nameV)]

  if(inherits(object, "glm") && type=="response") {
    mu<-object$family$linkinv(mu) 
  }
  #browser()
  if(interval!="none" || se.fit){
    if(!is.null(object$constr)){
      B <- do.call(blockdiag, list(diag(nLin), do.call(blockdiag, (object$constr$invA.RList)), diag(length(nameV))))
      V <- B %*% vcov(object) %*% t(B)
    } else {
      V <- vcov(object) 
      X <- X[,colnames(V)]
    }
    se <- sqrt(rowSums((X %*% V) * X))
    if(inherits(object, "glm")) {
      if(type=="response") se <- abs(object$family$mu.eta(eta))*se
      z<-abs(qnorm((1-level)/2)) 
      s2<-summary(object)$dispersion
    } else {
      z <- abs(qt((1-level)/2, df=object$df.residual))
      s2<- summary(object)$sigma^2
    }
    if(interval=="confidence"){
      mu<-cbind(fit=mu, lwr=mu-z*se, upr=mu+z*se)
    }
    if(interval=="prediction"){
      mu<-cbind(fit=mu, lwr=mu-z*sqrt(se^2+s2), upr=mu+z*sqrt(se^2+s2))
    }
  }
  if(se.fit) {
    mu <- list(fit=mu, se.fit=se, df= object$df.residual, residual.scale=sqrt(s2))
    if(!inherits(object, "glm")) mu$df<- object$df.residual
  }
  return(mu)
}
  
 
