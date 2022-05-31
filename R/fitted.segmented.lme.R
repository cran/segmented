fitted.segmented.lme<-function(object, level=1, sort=TRUE, ...){
  #fit: an object of class "segmented.lme" 
  #What about "fitted(oo$lme.fit.noG)" or "fitted(obj,level=1)+fit$Off"?
  #fitted(fitG,level=1)+fit$Off e' proprio uguale a fitted(fit.noG, level=1)  
  #comunque per level=0 (population parameter) l'identita' non vale, ed e' necessario fare i calcoli
  #   "manualmente"
  #obj<- object[[1]] #sarebbe fit$lme.fit
  levelC<- min(level,1) #valori >1 sono riportati ad 1
  levelC<-deparse(levelC)
  
  #browser()
  
  switch(levelC,
         "0"={
           #leftSlope<- if(object$namesGZ$nameZ %in% names(fixef(object[[2]]))) fixef(object[[2]])[object$namesGZ$nameZ] else 0
           #b0<-fixef(object[[2]])["(Intercept)"]
           #if(is.na(b0)) b0<-0
           r<-vector("list", length=length(names(object$psi.i)))
           mu0 <- fitted(object$lme.fit.noG, level=0)
           for(id in names(object$psi.i)){
             diffSlope<-object$fixed.eta.delta[paste(id)]
             Psi<- object$fixed.psi[paste(id)]
             x<- object$Z[names(object$Z)==id]
             #mu<-b0+leftSlope*x+diffSlope*pmax(x-Psi,0)
             psi.i <- object$psi.i[paste(id)]
             mu<-mu0[names(mu0)==paste(id)]
             mu <- mu  - diffSlope*pmax(x-psi.i,0) +diffSlope*pmax(x-Psi,0)
             r[[id]]<-mu
           }
           mu<-unlist(r)
           names(mu) <- names(object$Z)               
           #                mu<-fitted(obj,level=0) + fit$Off
           #                if("G0"%in%names(ranef(obj))){
           #                  ni<-tapply(obj$groups[,1], obj$groups[,1], length)
           #                  ki<-rep(ranef(obj)[["G0"]],ni)
           #                  mu<-mu + ki*obj$data[["G0"]]
           #                  }
         },
         "1"={ mu <- fitted(object$lme.fit.noG, level=level)
         #            "1"={ mu<-fitted(obj,level=1)+fit$Off #e' proprio uguale a fitted(fit[[2]], level=1)
         }
  ) #end_switch
  if(sort) mu<- mu[order(names(mu))]
  return(mu)
}

