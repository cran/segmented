getVarCov.segmented.lme<-function(obj, ..., type = c("random.effects", "conditional", 
                                                             "marginal")){
  type=match.arg(type)
  if(type=="marginal") warning("use type='marginal' at your own risk!")
  getVarCov(obj$lme.fit, ..., type=type)
}
