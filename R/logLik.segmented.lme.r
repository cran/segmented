logLik.segmented.lme<-function(object, ...){
  a<-logLik(object$lme.fit.noG)
  attr(a, "df") <-attr(logLik(object$lme.fit), "df")
  a
}
