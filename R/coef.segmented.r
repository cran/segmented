coef.segmented<-function(object, include.psi=FALSE, ...){
  b<- object$coefficients
  if(include.psi){
    psi<- object$psi[,"Est."]
    b[rownames(object$psi)]<-psi
  }
  b
}
