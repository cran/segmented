#new predict.segmented
predict.stepmented<-function(object, newdata, se.fit=FALSE, interval=c("none","confidence", "prediction"), 
                            type = c("link", "response"), na.action=na.omit, level=0.95, .coef=NULL, .vcov=NULL, 
                            apprx.fit=c("none","cdf"), apprx.se=c("cdf","none"), ...){
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
  
  dummy.matrix<-function(x.values, x.name, obj.seg, psi.est=TRUE, isV=FALSE, .coef=NULL, k=NULL){ 
    #given the segmented fit 'obj.seg' and a segmented variable x.name with corresponding values x.values,
    #this function simply returns a matrix with columns (x, (x-psi)_+, -b*I(x>psi))
    #if obj.seg does not include the coef for the linear "x", the returned matrix is  ((x-psi)_+, -b*I(x>psi)) 
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
    # nameU<-obj.seg$nameUV$U[f.U(obj.seg$nameUV$U,x.name)]
    # nameV<-obj.seg$nameUV$V[f.U(obj.seg$nameUV$V,x.name)] #grep(x.name, obj.seg$nameUV$V, value = TRUE)
    
    #browser()
    
    nameU<-object$nameUV$U
    nameV<-gsub("V","psi", object$nameUV$V)
    nameU<- nameU[f.U(nameU,x.name)]
    nameV<- nameV[f.U(nameV, x.name)] #grep(x.name, obj.seg$nameUV$V, value = TRUE)
    
    
    if(is.null(obj.seg$constr)){
      diffSlope<-estcoef[nameU]
    } else {
      diffSlope<-drop(obj.seg$constr$invA.RList[[match(x.name, obj.seg$nameUV$Z)]]%*%estcoef[nameU])[-1]
    }
    #browser()
    
    est.psi<-obj.seg$psi[nameV,"Est."]
    se.psi<-obj.seg$psi[nameV, "St.Err"]
    npsi <- length(est.psi)
    PSI <- matrix(est.psi, n, ncol = npsi, byrow=TRUE)
    minZ <- object$rangeZ[1, x.name]
    maxZ <- object$rangeZ[2, x.name]
    Z01<- (x.values-minZ)/(maxZ-minZ)
    PSI01 <- (PSI-minZ)/(maxZ-minZ)
    est.psi01<- (est.psi-minZ)/(maxZ-minZ)
    newd<-matrix(,length(x.values), length(nameU)+length(nameV))
    colnames(newd)<-c(sub("psi","U", nameV), nameV)
    #browser()
    for(j in 1:npsi){
      if(is.null(k)){
        idU<-match(nameU[j], nameU)
        snr.idU<-abs(estcoef[nameU][idU])/sigma(object)
        ss01=n^(-(.6 + .5*log(snr.idU)/sqrt(snr.idU) -abs(est.psi01[j]-.5)^(1/2)/length(object$residuals)^(1/2)))
        ss<- ss01*(maxZ-minZ)
        } else {
          ss=n^k
        }
      newd[ , nameU[idU]] <- pnorm((x.values-est.psi[j])/ss)
      newd[ , nameV[idU]] <-  -(estcoef[nameU][idU]/ss)*dnorm((x.values-est.psi[j])/ss)
      #newd<-cbind(x.values,dummy1,dummy2)
    }
    all.psi<-obj.seg$indexU[[x.name]]
    if(!is.null(all.psi) && length(all.psi)!=npsi){
      newZ<-matrix(x.values, length(x.values), npsi)
      nomi.psi.noti<-setdiff(names(all.psi),nameU)
      psi.noti<-setdiff(all.psi, est.psi)
      PSI.noti <- matrix(rep(psi.noti, rep(n, length(psi.noti))), ncol = length(psi.noti))
      nomi<-c(colnames(newd),nomi.psi.noti)
      newd<-cbind(newd, (newZ-PSI.noti)*(newZ>PSI.noti))
      colnames(newd)<-nomi
    }
    #browser()
    U<-sapply(est.psi, function(.x) 1*(x.values>.x))
    colnames(U) <- nameU
    newd<-list(U=U, newd=newd)
    #colnames(newd)[1]<-x.name
    #browser()
    return(newd)
  } #end dummy.matrix()
  estcoef <- if(is.null(.coef)) coef(object) else .coef
  if(is.null(names(estcoef))) stop("the coef estimates should be named")
  nLin<- length(setdiff(names(coef(object)), c(object$nameUV$U,object$nameUV$V)))
  nSeg<- length(object$nameUV$Z)
  type<-match.arg(type)
  interval<-match.arg(interval)
  apprx.fit <-match.arg(apprx.fit)
  apprx.se <-match.arg(apprx.se)
  if(inherits(object, "glm") && object$family$family!="gaussian" && interval=="prediction") 
    stop("prediction intervals are not allowed with non-gaussian glm")
  nameU<-object$nameUV$U
  nameV<-gsub("V","psi", object$nameUV$V)
  nameZ<-object$nameUV$Z
  
  #browser()
  
  if(missing(newdata)){
    X <- model.matrix.stepmented(object, type = apprx.se)
    X.noV <- model.matrix.stepmented(object, type = "no")
    colnomi.noV <-colnames(X.noV)
    if(apprx.fit=="cdf") X.noV[,nameU]<-X[,nameU]
    idNA<- rep(FALSE, nrow(X))
  } else {
    #browser()
    #nomiLin <- setdiff(all.vars(formula(object))[-1], c(object$nameUV$U,object$nameUV$V))
    nomiLin <- setdiff(all.vars(as.formula(paste("~",paste(formula(object))[3]))), c(nameU, nameV))
    if(any(is.na(match(nomiLin, names(newdata))))) stop(" 'newdata' should includes all variables")
    
    na.arg <- deparse(substitute(na.action))
    idNA<- !complete.cases(newdata)
    if(any(idNA)){
      newdata<-na.omit(newdata)
    }
    if(!na.arg%in%c("na.omit","na.pass")) stop("na.action should be 'na.omit' or 'na.pass'")
    
    
    n<-nrow(newdata)
    Ulist<-r<-NULL
    if(length(object$call$obj)>0){ #se l'ogg e' stato ottenuto da segmented.*
      
      # Fo<- formula(delete.response(terms(formula(eval(object$call$obj)))))
      # idSeg<- object$nameUV$Z %in% all.vars(Fo)
      # if(any(!idSeg)){
      #   Fo<- update.formula(Fo, as.formula(paste("~.+", paste(object$nameUV$Z[!idSeg], collapse="+"))))
      # }
      #nomiTerms, a differenza di nomiLin, include eventuali poly(w,2)
      nomiTerms<-setdiff(attr(terms(formula(object)),"term.labels"), c(nameU, nameV))
      idSeg<- object$nameUV$Z %in% nomiLin #potresti mettere anche "nomiTerms"
      if(any(!idSeg)){
        nomiTerms <- c(nomiTerms, object$nameUV$Z[!idSeg])
      }
      Fo<-as.formula(paste("~.+", paste(nomiTerms, collapse="+")))
      M<-model.matrix(Fo, data=newdata,
                      contrasts=object$contrasts, xlev = object$xlevels)
      } else { #se l'ogg e' stato ottenuto da stepreg
        
        #browser()
        
        Fo<-as.formula(object$nameUV$formulaSegAllTerms)
        if(any(all.vars(Fo)%in%names(object$xlevels))){
          M<-model.matrix(Fo, data=newdata, 
                        contrasts = object$contrasts, xlev=object$xlevels)
          } else {
            M<-model.matrix(Fo, data=newdata)
          }
        #nomiLin<- all.vars(object$formulaLin)[-1] #non funziona se la rispo e' cbind(y,n-y)
        nomiLin <- all.vars(as.formula(paste("~",paste(object$formulaLin)[3])))
        if(any(!nomiLin%in%all.vars(Fo))){
          #nomiLinOK<- nomiLin[!nomiLin%in%all.vars(Fo)]
          terminLin<-attr(terms(object$formulaLin),"term.labels")[!nomiLin%in%all.vars(Fo)]
          Fo <- as.formula(paste("~.-1+",paste(terminLin,collapse="+")))
          #Fo <- update.formula(Fo, as.formula(paste("~.+",paste(terminLin,collapse="+"))))
          M1<-model.matrix(Fo, data=newdata, 
                         contrasts = object$contrasts, xlev=object$xlevels)
          M<-cbind(M, M1) #[,nomiLinOK,drop=FALSE])
        }
      } #end se ogg e' stato ottenuto da segreg
    #browser()
    for(i in 1:length(nameZ)){
        x.values <- M[,nameZ[i]] 
        DM<- dummy.matrix(x.values, nameZ[i], object, k=list(...)$k)
        Ulist[[i]]<- DM$U
        r[[i]]<-DM$newd
    }
    #browser()
    X <-data.matrix(matrix(unlist(r), nrow=n, byrow = FALSE))
    colnames(X)<- unlist(sapply(r, colnames))
    X<-cbind(M,X)

    if("(Intercept)" %in% names(estcoef)) X<-cbind("(Intercept)"=1,X)
    #X<-X[,unique(colnames(X)),drop=FALSE]
    X<- X[, names(estcoef)]
    U<- data.matrix(matrix(unlist(Ulist), nrow=n, byrow = FALSE))
    colnames(U) <- nameU
    X.noV<-X
    X.noV[,nameU]<-U
    colnomi<- colnames(X)
    colnomi.noV <- setdiff(colnomi, nameV)
    X.noV <- X.noV[, colnomi.noV, drop=FALSE]
  } #end if non-missing(newdata)
  
  if(length(setdiff(colnames(X),names(estcoef)))>0) stop("error in the names (of the supplied newdata)")
  estcoef.noV<- estcoef[colnomi.noV]  

  #ignora eventuali altre variabili contenute in newdata
  #nomiOK<- intersect(names(estcoef.noV), colnames(X.noV))
  #X.noV<- X.noV[, nomiOK, drop=FALSE]
  #estcoef.noV<-estcoef.noV[nomiOK]
  
  mu <- eta<- drop(X.noV%*% estcoef.noV)
  
  if(!is.null(object$offset)) mu<- eta<- eta+ object$offset
  
  #ATTENZIONE c'e' il problema dell'appaiamento dei nomi!!!
  #il problema e' che estcoef non ha sempre nomi!! 
  
  if(inherits(object, "glm") && type=="response") {
    mu<-object$family$linkinv(mu) 
  }
  
  if(!se.fit && interval=="none"){
    if(any(idNA) && na.arg=="na.pass"){
      mu0<-mu
      mu<- rep(NA, length(idNA))
      mu[!idNA]<-mu0
    }
    return(mu)
  } else { # se if(interval!="none" || se.fit)
    #browser()
    V <- if(is.null(.vcov)) vcov.stepmented(object, type=apprx.se, ...) else .vcov
    
    if(!is.null(object$constr)){
      B=if(nLin>0) append(list(diag(nLin)), object$constr$invA.RList, 1) else object$constr$invA.RList
      B=append(B, list(diag(length(nameV))), 2)
      B= do.call(blockdiag, B)
      V <- B %*% V %*% t(B)
    } else {
      X <- X[,colnames(V)] #semplicemente elimina e riordina le colonne di X
    }
    se <- sqrt(rowSums((X %*% V) * X))
    if(inherits(object, "glm")) {
      if(type=="response") se <- abs(object$family$mu.eta(eta))*se
      z<-abs(qnorm((1-level)/2)) 
      s2<- sigma(object)^2 #summary(object)$dispersion
    } else {
      z <- abs(qt((1-level)/2, df=object$df.residual))
      s2<- sigma(object)^2 #summary(object)$sigma^2
    }
    
    if(any(idNA) && na.arg=="na.pass"){
      mu0<-mu
      se0<-se
      mu<-se<- rep(NA, length(idNA))
      mu[!idNA]<-mu0
      se[!idNA]<-se0
    }
    
    
    if(interval=="confidence"){
      mu<-cbind(fit=mu, lwr=mu-z*se, upr=mu+z*se)
    }
    if(interval=="prediction"){
      mu<-cbind(fit=mu, lwr=mu-z*sqrt(se^2+s2), upr=mu+z*sqrt(se^2+s2))
    }
  
  if(se.fit) {
    mu <- list(fit=mu, se.fit=se, df= object$df.residual, residual.scale=sqrt(s2))
    if(!inherits(object, "glm")) mu$df<- object$df.residual
  }
  return(mu)
}
}

