fixef.segmented.lme<-fixed.effects.segmented.lme<-function(object,...) {
  b.all<-object$lme.fit$coefficients$fixed
  b.noG<-object$lme.fit.noG$coefficients$fixed
  b.all[names(b.noG)]<- b.noG
  b.all
}
