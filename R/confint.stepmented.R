`confint.stepmented` <- function(object, parm, level=0.95, method=c("delta", "score", "gradient"), round=TRUE, cheb=FALSE,
                                #var.diff=FALSE
                                digits=max(4, getOption("digits") - 1), .coef=NULL, .vcov=NULL, ...){
  method<-match.arg(method)
  if(method!="delta") stop("Only delta allowed")
  
  if(missing(parm)) {
    nomeZ<- object$nameUV$Z
  } else {
    if(is.numeric(parm)) parm<-object$nameUV$Z[parm]
    if(! all(parm %in% object$nameUV$Z)) stop("invalid 'parm' name", call.=FALSE)
    nomeZ<-parm
  }

  if(length(nomeZ)>1) {
    warning("There are multiple segmented terms. The first is taken", call.=FALSE, immediate. = TRUE)  
    nomeZ<-nomeZ[1]
  }

  nomiZ<- object$nameUV$Z
  nomiV<- object$nameUV$V
  nomiU<- object$nameUV$U
  nomiPsi <- gsub("V","psi", nomiV)
  
  Cov<-vcov.stepmented(object,...)
  id <- match(nomiPsi, names(coef(object)))
  vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
  
  psi<-object$psi[nomiPsi,"Est."]
  se<- sqrt(vv) #[nomiPsi]
  npsi<-length(psi)
  
  if(cheb){
    z<-1/sqrt(1-level)
  } else {
    z<-if("lm"%in%class(object)) abs(qt((1-level)/2,df=object$df.residual)) else abs(qnorm((1-level)/2))
  }
  #browser()
  #z=abs(qnorm((1-level)/2))
  Z<-object$Z
  Z0<-apply(Z,2,sort)
  
  #browser()
  
  inf<-pmax(psi -z*se, object$rangeZ[1,])
  sup<-pmin(psi +z*se, object$rangeZ[2,])
  
  #ripeti i nomi delle variabili stepmented tante volte quanti sono i psi..
  nomiZripetuti<- sub("\\.", "", sub("psi[1-9].","", nomiPsi))
  if(round){
    inf.rounded<-sapply(1:npsi, function(j) Z0[sum(Z0[, nomiZripetuti[j]]<inf[j])+c(0,1),nomiZripetuti[j]])
    sup.rounded<-sapply(1:npsi, function(j) Z0[sum(Z0[, nomiZripetuti[j]]<sup[j])+c(0,1),nomiZripetuti[j]])
    r<-cbind(object$psi.rounded[1,], inf.rounded[1,], sup.rounded[1,])
  } else {
    r<-cbind(psi, inf, sup)
  }
  
  rownames(r)<-nomiPsi
  colnames(r)<- c("Est.",paste("CI","(",level*100,"%",")",c(".low",".up"),sep=""))
  r<-r[endsWith(nomiPsi, paste(".", nomeZ ,sep="")),]
  r
}

