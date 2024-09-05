#new predict.segmented
predict.segmented<-function(object, newdata, se.fit=FALSE, interval=c("none","confidence", "prediction"), 
            type = c("link", "response"), na.action=na.omit,# "terms"),
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
    #le seguenti righe selezionavano (ERRONEAMENTE) sia "U1.x" sia "U1.neg.x" (se "x" e "neg.x" erano segmented covariates)
    #nameU<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$U, value = TRUE)
    #nameV<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$V, value = TRUE)
    nameU<-obj.seg$nameUV$U[f.U(obj.seg$nameUV$U,x.name)]
    nameV<-obj.seg$nameUV$V[f.U(obj.seg$nameUV$V,x.name)] #grep(x.name, obj.seg$nameUV$V, value = TRUE)
    
    
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
  if(is.null(names(estcoef))) stop("the coef estimates should be named")
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
  
  #browser()
  
  if(missing(newdata)){
    X <- model.matrix(object)
    idNA<- rep(FALSE, nrow(X))
  } else {
    #browser()
    #nomiLin <- setdiff(all.vars(formula(object))[-1], c(object$nameUV$U,object$nameUV$V))
    nomiLin <- setdiff(all.vars(as.formula(paste("~",paste(formula(object))[3]))), c(object$nameUV$U,object$nameUV$V))
    if(any(is.na(match(nomiLin, names(newdata))))) stop(" 'newdata' should includes all variables")
    #devi trasformare la variabili segmented attraverso dummy.matrix()
    
    na.arg <- deparse(substitute(na.action))
    idNA<- !complete.cases(newdata)
    if(any(idNA)){
      newdata<-na.omit(newdata)
    }
    if(!na.arg%in%c("na.omit","na.pass")) stop("na.action should be 'na.omit' or 'na.pass'")
    
    n<-nrow(newdata)
    r<-NULL
    if(length(object$call$obj)>0){ #se l'ogg e' stato ottenuto da segmented.*
      # Fo<- formula(delete.response(terms(formula(eval(object$call$obj)))))
      # idSeg<- object$nameUV$Z %in% all.vars(Fo)
      # if(any(!idSeg)){
      #   Fo<- update.formula(Fo, as.formula(paste("~.+", paste(object$nameUV$Z[!idSeg], collapse="+"))))
      # }
      
      #nomiTerms, a differenza di nomiLin, include eventuali poly(w,2)
      nomiTerms<-setdiff(attr(terms(formula(object)),"term.labels"),c(object$nameUV$U,object$nameUV$V))
      idSeg<- object$nameUV$Z %in% nomiLin #potresti mettere anche "nomiTerms"
      if(any(!idSeg)){
        nomiTerms <- c(nomiTerms, object$nameUV$Z[!idSeg])
      }
      
      Fo<-as.formula(paste("~.+", paste(nomiTerms, collapse="+")))
      
      M<-model.matrix(Fo, data=newdata,
              contrasts=object$contrasts, xlev = object$xlevels)
      
    } else { #se l'ogg e' stato ottenuto da segreg
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

    }
    
    for(i in 1:length(nameZ)){
      x.values <- M[,nameZ[i]] 
      DM<-dummy.matrix(x.values, nameZ[i], object)
      r[[i]]<-DM
    }
    
    #browser()
    X <-data.matrix(matrix(unlist(r), nrow=n, byrow = FALSE))
    colnames(X)<- unlist(sapply(r, colnames))
    X<-cbind(M,X)
    X<-X[,unique(colnames(X)),drop=FALSE]
    if("(Intercept)" %in% names(estcoef)) X<-cbind("(Intercept)"=1,X)
  }


  if(!is.null(object$constr)){
    for(i in 1:length(nameZ)){
        nomeU.i<-grep(object$nameUV$Z[i], object$nameUV$U, value=TRUE)
        idU.i <- match(nomeU.i, names(estcoef))
        coef.new<-drop(object$constr$invA.RList[[i]]%*%estcoef[nomeU.i])
        names(coef.new)<-c(object$nameUV$Z[i], 
                           paste("U",1:(length(coef.new)-1),".",object$nameUV$Z[i],sep="" ))
        estcoef<-append(estcoef[-idU.i], coef.new, after=idU.i[1]-1)
    }
  }
  X<-X[,names(estcoef),drop=FALSE]

  if(length(setdiff(colnames(X),names(estcoef)))>0) stop("error in the names (of the supplied newdata)")
  
  #browser()
  colnomi<- colnames(X)
  
  colnomi.noV <- setdiff(colnomi, nameV)
  X.noV <- X[, colnomi.noV, drop=FALSE]
  estcoef.noV<- estcoef[colnomi.noV]

  #ignora eventuali altre variabili contenute in newdata
  #nomiOK<- intersect(names(estcoef.noV), colnames(X.noV))
  #X.noV<- X.noV[, nomiOK, drop=FALSE]
  #estcoef.noV<-estcoef.noV[nomiOK]
  
  mu <- eta<- drop(X.noV%*% estcoef.noV)
  
  if(!is.null(object$offset)) mu<- eta<- eta+ object$offset

  X <- X[,c(colnomi.noV, nameV),drop=FALSE]

  if(inherits(object, "glm") && type=="response") {
    mu<-object$family$linkinv(mu) 
  }
  #browser()
  if(interval!="none" || se.fit){
    V <- vcov(object) 
    if(!is.null(object$constr)){
      B=if(nLin>0) append(list(diag(nLin)), object$constr$invA.RList, 1) else object$constr$invA.RList
      B=append(B, list(diag(length(nameV))), 2)
      B= do.call(blockdiag, B)
      V <- B %*% V %*% t(B)
    } else {
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
  } else {
    if(any(idNA)&& na.arg=="na.pass"){
      mu0<-mu
      mu<- rep(NA, length(idNA))
      mu[!idNA]<-mu0
    }
  }
  if(se.fit) {
    mu <- list(fit=mu, se.fit=se, df= object$df.residual, residual.scale=sqrt(s2))
    if(!inherits(object, "glm")) mu$df<- object$df.residual
  }
  return(mu)
}
  
 
