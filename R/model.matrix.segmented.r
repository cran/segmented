model.matrix.segmented<-function(object,...){
  #if(!inherits(object, "segmented")) stop("A 'segmented' fit is requested")
  if(inherits(object, "lm")) {
    X<- qr.X(object$qr, ...)
    if(inherits(object, "glm")) {
      #W<-chol(diag(object$weights))
      #X <- X/diag(W)
      X<- X/sqrt(object$weights)
    }
  } else {
    class(object)<-class(object)[-1]
    X<-try(model.matrix(object,...), silent=TRUE)
    if(!is.matrix(X)) X<- model.matrix(object, data=model.frame(object))
  }
  X
}
