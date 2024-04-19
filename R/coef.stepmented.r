coef.stepmented <-function(object, include.psi=TRUE, ...){
  #browser()
  b<- object$coefficients #solo coeffs lineari (senza psi)
  if(!all(match( colnames(object$psi.rounded), names(b), 0)>0)) {
    psi<- object$psi.rounded[1,]
    names(psi)<-colnames(object$psi.rounded)
    b <- c(b, psi)
  }
  if(!include.psi) {
    b[match( colnames(object$psi.rounded), names(b), 0)] <-0
  }
  b
}
