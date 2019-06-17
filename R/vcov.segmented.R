#vc<-function(obj){
#  invXtX<-chol2inv(qr.R(obj$qr)) #(XtX)^{-1}
#  V<-vcov.segmented(obj,is=TRUE) 
#  s2<- if(inherits(obj, "glm")) summary.glm(obj)$dispersion else summary.lm(obj)$sigma^2
#  s2*V%*% invXtX %*% V
#}

vcov.segmented<-function(object, var.diff=FALSE, is=FALSE, ...){
  if(is && inherits(object, "Arima")) {
    warning("is=TRUE ignored with Arima fits", call.=FALSE)
    is<-FALSE
  }
  if(is){
    if(var.diff) warning("option 'var.diff=TRUE' ignored with 'is=TRUE' ", call.=FALSE)
    X<-model.matrix(object) #qr.X(object$qr) piu efficiente?
    nomiZ<- object$nameUV$Z
    nomiV<- object$nameUV$V
    nomiU<- object$nameUV$U
    for(i in 1:length(nomiV)){
      nomeU<-nomiU[i]
      nomeV<-nomiV[i]
      nomepsi<-strsplit(nomeV,"\\.")[[1]][1] #solo "psi1" o "psi2",.. e' meglio estrarre il "psi1" perche' il nome della variabile puo' contenere un punto.. 
      nomeZ<-gsub(paste(nomepsi,".",sep=""),"",nomeV) #estrae il nome della variabile..
      Z<-X[,nomeZ]
      est.psi<- object$psi[nomeV,"Est."]
      se.psi<- object$psi[nomeV,"St.Err"]
      X[,nomeV]<- (-object$coefficients[nomeU])*pnorm((Z-est.psi)/se.psi)
      }
    s2<- if(inherits(object, "glm")) summary.glm(object)$dispersion else summary.lm(object)$sigma^2
    w<-object$weights
    if(is.null(w)) w<-1
    v<-s2*solve(crossprod(X*sqrt(w)))
    return(v)
      } else {
    if(inherits(object, "Arima")){ 
      v<-object$var.coef
      return(v)
    }
    if(inherits(object, "glm")){
        if(var.diff) warning("option 'var.diff=TRUE' ignored with 'glm' objects", call.=FALSE)
        so <- summary.glm(object, correlation = FALSE, ...)
        v<-so$dispersion * so$cov.unscaled
        return(v)
    } 
    if(inherits(object, "lm")){
        if(var.diff){
            if(length(object$nameUV$Z)>1) {
                var.diff<-FALSE
                warning("var.diff set to FALSE with multiple segmented variables", call.=FALSE)
            }
          v<-summary.segmented(object, var.diff=TRUE, correlation = FALSE, ...)$cov.var.diff
        } else {
          so<-summary.segmented(object, var.diff=FALSE, correlation = FALSE, ...)
          v<-so$sigma^2 * so$cov.unscaled #object$cov.unscaled.is
          }
        return(v)
    } else { #in tutti gli altri casi..
      if(class(object)[1]=="segmented") class(object)<-class(object)[-1]
      v<-vcov(object)
      #paste("vcov.",class(object),sep="")
      return(v)
    }
  } #end else is
} #end fn



