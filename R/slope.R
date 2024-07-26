`slope` <- function(ogg, parm, conf.level=0.95, rev.sgn=FALSE, APC=FALSE, .vcov=NULL, .coef=NULL, use.t=NULL,
                    by=NULL, interc=TRUE, ..., digits = max(4, getOption("digits") - 2)){

#   se e' un "newsegmented"
#    if(!is.null(ogg$R.slope)) {
#             covv<-old.coef.var(ogg)
#             ogg$coefficients<- covv$b
#             covv<-  covv$cov
#             ogg$psi<-old.psi(ogg)
#             ogg$nameUV<-old.nomi(ogg)
#             } else {
#             covv<-try(vcov(ogg,...), silent=TRUE)
#             }
  slopeM<-function(obj, by=NULL, conf.level=0.95, vcov.=NULL, ... ){
    #=========>> da provare con by con piu' termini e se leftSlope=0 
    #obj: the segmented.lme object
    #by: a named list indicating covariate names and corresponding values affecting the fitted segmented relationship.
    #    Example: by=list(group="2",z2=.2) 
    #conf.level: for pointwise CI..
    #withI: if TRUE, the fitted lines are plotted with intercept (if included in the model)
    #vcov.: the fixed effect cov matrix. If NULL is computed by vcov.segmented.lme
    #drop.var: possible coefficient names to be removed before computing the segmented relationship (E.g. the group-specific intercept..)
    #obj[[1]] -> obj$lme.fit
    #obj[[2]] -> obj$lme.fit.noG
    #------------------
    V<- if(is.null(vcov.)) vcov.segmented.lme(obj) else vcov. #object$lme.fit$varFix 
    if(!is.null(by) && !is.list(by)) stop("if provided, 'by' should be a (named) list of scalars")
    Z<-obj$Z 
    nomeZ<-obj$namesGZ$nameZ
    beta.noG<- fixef(obj$lme.fit.noG) 
    beta.all<-fixef(obj$lme.fit)
    #browser()
    
    beta.G<-beta.all[setdiff(names(fixef(obj$lme.fit)), names(beta.noG))]
    nomiCoef<-names(beta.noG)
    if(!is.null(by)) {
      a<-by
      #isZero<-sapply(a, function(x) x==0)
      if(!all(sapply(a, length)==1)) stop("vectors in 'by' are not allowed")
      nomiOK<-const<-idList<-vector("list", length(a))
      values<-vector(,length(a))
      
      for(i in 1:length(a)) {
        nomiOK[i]<-nomeOK <- if(is.character(a[[i]])) paste(names(a[i]),a[[i]], sep="") else names(a[i])
        
        #replace 0
        if(a[[i]]==0) a[[i]]<- 1e-16
        
        #per la left slope
        bLeftSlope<-c(beta.noG[paste(obj$namesGZ$nameZ,":",nomeOK, sep="")],
                      beta.noG[paste(nomeOK, ":", obj$namesGZ$nameZ, sep="")])
        bLeftSlope<-bLeftSlope[!is.na(bLeftSlope)]
        if(!is.character(a[[i]])) bLeftSlope<-bLeftSlope*a[[i]]
        if(length(bLeftSlope)<=0) bLeftSlope<-NA 
        
        #per la slope-diff
        bU<-beta.noG[paste("U", nomeOK, sep=".")]
        if(!is.character(a[[i]])) bU <-bU*a[[i]]
        
        const[[i]]<-c(bLeftSlope, bU)
        const[[i]]<- ifelse(is.na(const[[i]]),0,const[[i]])
        
        idList[[i]]<-names(c(bLeftSlope, bU))
        values[i]<-ifelse(is.character(a[[i]]),1,a[[i]])
      }
      #browser()          
      const<-matrix(unlist(const),2, byrow=FALSE)
      colnames(const)<-names(by)
      nomiNOdiff <- names(which(colSums(const)==0))
      if(length(nomiNOdiff)>0) warning("The value of", paste(" '", paste(nomiNOdiff, collapse=" and "),"' ",sep=""), 
                                       "supplied in 'by' does not modify the baseline slopes estimates", call. = FALSE)
      
      nomiCoef<- c(nomeZ, "U", unlist(idList))
      ##########################################    
    } else { #se 'by' e' NULL
      const<-matrix(0,2,1) 
      nomiCoef<- c(nomeZ, "U")
      values<-c(1,1)
    }  
    ##########################################    
    #browser()  
    final.names<-setdiff(nomiCoef, c("G0",obj$namesGZ$nomiG,""))
    final.names<-final.names[!is.na(final.names)]
    
    #prepara la matrice del disegno..
    X<-matrix(1, 2, 2) #SOLO un breakpoint, quindi 2 segmenti..
    X[row(X)<col(X)]<-0
    Ident<-diag(ncol(X))
    M<-vector("list", length=ncol(const))
    for(j in 1:ncol(const)) M[[j]]<-values[j]*Ident[, which(const[,j]!=0), drop=FALSE]
    M<-cbind(Ident, do.call("cbind", M))
    XX<- X%*%M
    r<- drop(XX %*% beta.noG[final.names])
    V<-V[final.names,final.names]
    SE.fit<-sqrt(rowSums((XX %*% V) * XX)) #sqrt(diag(X%*%Var%*%t(X)))
    zalpha<- -qnorm((1-conf.level)/2)
    ci<-cbind(r-zalpha*SE.fit, r+zalpha*SE.fit)
    #colnames(ci)<-paste("CI", "(", conf.level * 100,"%", ")", c(".low", ".up"), sep = "")
    colnames(ci)<-paste(" ", conf.level, c(".low", ".up"), sep = "")
    ris<-cbind("Est."=r, "St.Err"=SE.fit, "t value"=r/SE.fit, ci)
    rownames(ris)<-c("leftSlope", "rightSlope")
    #ris
    #  rev.sgn=FALSE
    #  if(rev.sgn){
    #    ris<-cbind(-ris[nrow(ris):1,1],ris[nrow(ris):1,2],-ris[nrow(ris):1,3],
    #               -ris[nrow(ris):1,5],-ris[nrow(ris):1,4])}
    ris
  }
  
  #se metti la linea sotto, poi NON funziona con oggetti che provengono da segmented.default 
  #if(!inherits(ogg, "stepmented") || !inherits(ogg, "segmented") || !inherits(ogg, "segmented.lme")) stop("slope() works only for..")
  
  if(inherits(ogg, "stepmented")){
    covv <- if(is.null(.vcov)) vcov(ogg, ...) else .vcov 
    estcoef<- if(is.null(.coef)) coef(ogg) else .coef
    if(length(estcoef)==0) stop("No coefficient in the object fit?")
    if(!all(dim(covv)==c(length(estcoef), length(estcoef)))) stop("dimension of cov matrix and estimated coeffs do not match", call. = FALSE)
    nomepsi<-rownames(ogg$psi) #OK
    nomeU<-ogg$nameUV$U
    nomeZ<-ogg$nameUV$Z
    if(missing(parm)) {
      nomeZ<- ogg$nameUV$Z
    } else {
      if(! all(parm %in% ogg$nameUV$Z)) {
        stop("invalid parm") } else {nomeZ<-parm}
    }
    nomi<-names(estcoef)
    index<-vector(mode = "list", length = length(nomeZ))
    
    for(i in 1:length(nomeZ)) {
      idU <- ogg$nameUV$U[endsWith(ogg$nameUV$U, paste(".", nomeZ[i], sep = ""))]
      if (interc && "(Intercept)" %in% nomi) idU <- c("(Intercept)", idU)
      index[[i]]<- idU
    }
    
    if(!is.null(use.t)){
      k0<- if(use.t) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
    } else {
      k0<-if(inherits(ogg, "lm", which=TRUE)==2) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
    }
    #k<-if("lm"%in%class(ogg)) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
    Ris<-list()   
    #digits <- max(3, getOption("digits") - 3)
    for(i in 1:length(index)){
      #ind<-as.numeric(na.omit(unlist(index[[i]])))
      ind<- match(na.omit(unlist(index[[i]])),nomi)
      M<-matrix(1,length(ind),length(ind))
      M[row(M)<col(M)]<-0
      cof<-estcoef[ind]
      cof.out<-M%*%cof 
      
      if(!inherits(covv, "try-error")){
        cov.ok<-covv[ind,ind]
        cov.out<-M%*%cov.ok%*%t(M)
        se.out<-sqrt(diag(cov.out))
        k<-k0*se.out
        ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
      } else {
        ris<-cbind(cof.out, NA, NA, NA,NA)                
      }
      cin<-paste("CI","(",conf.level*100,"%",")",c(".l",".u"),sep="")
      #se la left slope e' nulla... per stepmented NON vale!
      #if(!nomeZ[i]%in%nomi){            
      #  ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)
      #}
      if(!interc || !("(Intercept)"%in%nomi) ){
        ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)
      }
      nomeT<-if("lm"%in%class(ogg)) "t value" else "z value"
      dimnames(ris)<-list(paste("level", 1:nrow(ris), sep=""),c("Est.","St.Err.",nomeT,cin[1],cin[2]))
      Ris[[nomeZ[i]]]<- signif(ris,digits)
      
    } #end loop i
    return(Ris)
  } else {
    if(length(ogg)==2 && all(sapply(ogg, function(.x)inherits(.x, "segmented") ))){
      obj1<-ogg[[1]]
      obj2<-ogg[[2]]
      os1 <- slope(obj1,conf.level=conf.level,...)[[1]]
      os2 <- slope(obj2,conf.level=conf.level,...)[[1]]
      K<-min(length(os1[,1]),length(os2[,1]))
      estdiff<-zvalue<-pvalue<-rep(NA,K)
      for(id in 1:K){
        est1<- os1[id,1]
        se1 <- os1[id,2]
        est2<- os2[id,1]
        se2<- os2[id,1]
        estdiff[id]<-abs(est1-est2)
        zvalue[id] <-(est1-est2)/sqrt(se1^2+se2^2)
        pvalue[id] <-2*(1-pnorm(abs(est1-est2)))
      }
      r<-cbind("|est diff|"=estdiff, "p-value"=signif(pvalue,3))
      rownames(r)<-paste("slope",1:K,sep="")
      #if(msg) 
      cat("Comparing 'matched' slopes for two segmented fits\n")
      r
    } else {
    #if(length(ogg)>2 && (inherits(ogg, "segmented") || inherits(ogg, "segmented.lme"))){
      if(class(ogg)[1]=="segmented.lme"){
            a<-slopeM(ogg, conf.level=conf.level, vcov.=.vcov, by=by, ...)
            return(a)
            
          } else {
            covv <- if(is.null(.vcov)) vcov(ogg, ...) else .vcov 
            if(is.null(.coef)) {
              estcoef<- coef(ogg) 
              if(is.null(estcoef)) estcoef <- ogg$coef
              if(is.null(estcoef)) stop("No coeffs in the fit? Please use '.coef'")
            } else {
              estcoef<- .coef
            }
            #browser()
            if(length(estcoef)==0) stop("No coefficient in the object fit?")
    
            if(!all(dim(covv)==c(length(estcoef), length(estcoef)))) stop("dimension of cov matrix and estimated coeffs do not match", call. = FALSE)
            
            nomepsi<-rownames(ogg$psi) #OK
            nomeU<-ogg$nameUV$U
            nomeZ<-ogg$nameUV$Z
            if(missing(parm)) {
              nomeZ<- ogg$nameUV$Z
              if(length(rev.sgn)==1) rev.sgn<-rep(rev.sgn,length(nomeZ))
              } else {
                if(! all(parm %in% ogg$nameUV$Z)) {
                  stop("invalid parm") } else {nomeZ<-parm}
                }
            if(length(rev.sgn)!=length(nomeZ)) rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
            nomi<-names(estcoef)
            index<-vector(mode = "list", length = length(nomeZ))
            for(i in 1:length(nomeZ)) {
              index[[i]]<- match(c(nomeZ[i],ogg$nameUV$U[grep(nomeZ[i], ogg$nameUV$U)]), names(estcoef),0)
              n.psi.est.i <- length(ogg$nameUV$V[grep(nomeZ[i], ogg$nameUV$V)])
              if(length(ogg$indexU[[nomeZ[i]]])!=n.psi.est.i){ #se ci sono anche psi fissi
                if(!is.null(ogg$constr)){
                  stop("slope() does not work with constraints and fixed psi")
                } else {
                  index[[i]]<-match(c(nomeZ[i], names(ogg$indexU[[nomeZ[i]]])),  nomi, 0)
                }
              }
                
            }
            if(!is.null(use.t)){
              k0<- if(use.t) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
            } else {
              k0<-if(inherits(ogg, "lm", which=TRUE)==2) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
            }
            #k<-if("lm"%in%class(ogg)) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
            Ris<-list()   
            #digits <- max(3, getOption("digits") - 3)
            rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
            
            #browser()
            
            for(i in 1:length(index)){
                ind<- index[[i]] #as.numeric(unlist(index[[i]]))
                ind<-ind[ind!=0]
                cof<-estcoef[ind]
                if(is.null(ogg$constr)){
                  M<-matrix(1,length(ind),length(ind))
                  M[row(M)<col(M)]<-0
                  
                  cof.out<-M%*%cof
                } else {
                  M<-ogg$constr$RList[[match(nomeZ[i], ogg$nameUV$Z,0)]]
                  cof.out <- M%*%cof
                }
    
                if(!inherits(covv, "try-error")){
                    cov.ok<-covv[ind,ind]
                    cov.out<-M%*%cov.ok%*%t(M)
                    se.out<-sqrt(diag(cov.out))
                    k<-k0*se.out
                    ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
                    } else {
                    ris<-cbind(cof.out, NA, NA, NA,NA)                
                    }
                    cin<-paste("CI","(",conf.level*100,"%",")",c(".l",".u"),sep="")
                #se la left slope e' nulla....
                #if(identical(length(ind),length(grep(paste("\\.",nomeZ[i],"$",sep=""), nomeU)))){
                if(is.null(ogg$constr) && !nomeZ[i]%in%nomi){
                  ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)
                }
                    
                if(rev.sgn[i]){
                    ris<-cbind(-ris[nrow(ris):1,1],ris[nrow(ris):1,2],-ris[nrow(ris):1,3],
                          -ris[nrow(ris):1,5],-ris[nrow(ris):1,4])
                          }
                nomeT<-if("lm"%in%class(ogg)) "t value" else "z value"
                dimnames(ris)<-list(paste("slope", 1:nrow(ris), sep=""),c("Est.","St.Err.",nomeT,cin[1],cin[2]))
                if(APC) ris<-100*(exp(ris[,c(1,4,5)])-1)
                Ris[[nomeZ[i]]]<- signif(ris,digits)
                
                    } #end loop i
                #if(!missing(parm)){
                #    if(!all(parm %in% ogg$nameUV[[3]])) stop("invalid parm") else Ris<-Ris[parm]}
                return(Ris)
          }
    }
  }
}


  
  # else {
  #   obj1<-ogg[[1]]
  #   obj2<-ogg[[2]]
  #   os1 <- slope(obj1,conf.level=conf.level,...)[[1]]
  #   os2 <- slope(obj2,conf.level=conf.level,...)[[1]]
  #   K<-min(length(os1[,1]),length(os2[,1]))
  #   estdiff<-zvalue<-pvalue<-rep(NA,K)
  #   for(id in 1:K){
  #     est1<- os1[id,1]
  #     se1 <- os1[id,2]
  #     est2<- os2[id,1]
  #     se2<- os2[id,1]
  #     estdiff[id]<-abs(est1-est2)
  #     zvalue[id] <-(est1-est2)/sqrt(se1^2+se2^2)
  #     pvalue[id] <-2*(1-pnorm(abs(est1-est2)))
  #   }
  #   r<-cbind("|est diff|"=estdiff, "p-value"=signif(pvalue,3))
  #   rownames(r)<-paste("slope",1:K,sep="")
  #   #if(msg) 
  #   cat("Comparing 'matched' slopes for two segmented fits\n")
  #   r
  # }
  


