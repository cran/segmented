breakpoint<-function(ogg, by=NULL, conf.level=0.95, vcov.=NULL, use.link=TRUE, ... ){
  #=========>> da provare con by con piu' termini e se leftSlope=0 
  #ogg: the segmented.lme oggect
  #by: a named list indicating covariate names and corresponding values affecting the fitted segmented relationship.
  #    Example: by=list(group="2",z2=.2) 
  #conf.level: for pointwise CI..
  #vcov.: the fixed effect cov matrix. If NULL is computed by vcov.segmented.lme
  #drop.var: possible coefficient names to be removed before computing the segmented relationship (E.g. the group-specific intercept..)
  #ogg[[1]] -> ogg$lme.fit
  #ogg[[2]] -> ogg$lme.fit.noG
  #------------------
  #trasformare evntualmente se link logit
  #calcolare pure un mu(psi)?
  #-------------------
  #if(level>0) stop("only level=0 allowed")
  
  #questi sono i codici presi da plot.segmented.lme()
  # Z<-obj$Z 
  # beta.noG<- fixef(obj$lme.fit.noG) 
  # beta.all<-fixef(obj$lme.fit)
  # beta.G<-beta.all[setdiff(names(fixef(obj$lme.fit)), names(beta.noG))]
  # nomiCoef<-names(beta.noG)
  # nomeZ<-obj$namesGZ$nameZ
  # nomiCoef <- if(nomeZ%in%nomiCoef) c(nomeZ, "U", "G0") else c("U", "G0")
  # if("(Intercept)"%in%names(beta.noG)) nomiCoef <-c("(Intercept)", nomiCoef)
  # 
  # 
  # #browser()
  # 
  # if(!is.null(by)) {
  #   a<-by
  #   #isZero<-sapply(a, function(x) x==0)
  #   if(!all(sapply(a, length)==1)) stop("vectors in 'by' are not allowed")
  #   nomiOK<-const<-idList<-vector("list", length(a))
  #   values<-vector(,length(a))
  #   for(i in 1:length(a)) {
  #     nomiOK[i]<-nomeOK <- if(is.character(a[[i]])) paste(names(a[i]),a[[i]], sep="") else names(a[i])
  #     
  #     #replace 0
  #     if(a[[i]]==0) a[[i]]<- 1e-16
  #     
  #     #per l'intercetta
  #     bInterc<- beta.noG[nomeOK]
  #     if(!is.character(a[[i]])) bInterc<-bInterc*a[[i]] 
  #     
  #     #per la left slope
  #     bLeftSlope<-c(beta.noG[paste(obj$namesGZ$nameZ,":",nomeOK, sep="")],
  #                   beta.noG[paste(nomeOK, ":", obj$namesGZ$nameZ, sep="")])
  #     bLeftSlope<-bLeftSlope[!is.na(bLeftSlope)]
  #     if(!is.character(a[[i]])) bLeftSlope<-bLeftSlope*a[[i]]
  #     if(length(bLeftSlope)<=0) bLeftSlope<-NA 
  #     
  #     #per la slope-diff
  #     bU<-beta.noG[paste("U", nomeOK, sep=".")]
  #     if(!is.character(a[[i]])) bU <-bU*a[[i]]
  #     
  #     #per il changepoint
  #     bG<-beta.G[paste("G", nomeOK, sep=".")]
  #     if(!is.character(a[[i]])) bG <-bG*a[[i]]
  #     
  #     const[[i]]<-c(bInterc, bLeftSlope, bU, bG)
  #     const[[i]]<- ifelse(is.na(const[[i]]),0,const[[i]])
  #     
  #     idList[[i]]<-names(c(bInterc, bLeftSlope, bU, bG))
  #     values[i]<-ifelse(is.character(a[[i]]),1,a[[i]])
  #   }
  #   #browser()          
  #   const<-matrix(unlist(const),4, byrow=FALSE)
  #   colnames(const)<-names(by)
  #   nomiNOdiff <- names(which(colSums(const)==0))
  #   if(length(nomiNOdiff)>0) warning("The", paste(" '", nomiNOdiff,"' ",sep=""), "value supplied in 'by' does not modify the baseline line", call. = FALSE)
  #   
  #   nomiCoef<- c(nomiCoef, unlist(idList))
  #   ##########################################
  # } else { #se 'by' e' NULL
  #   const<-matrix(0, length(nomiCoef),1)
  #   values<-rep(1, length(nomiCoef))
  # }
  # ##########################################    
  # #browser()  
  # 
  # #prepara la matrice del disegno..
  # est.psi.fixed <-fixef(obj$lme.fit)["G0"]+ sum(const[4,])#sum(const[length(nomiCoef),])
  
  
  
  
  
  V<- if(is.null(vcov.)) vcov.segmented.lme(ogg,...) else vcov. #oggect$lme.fit$varFix 
  if(!is.null(by) && !is.list(by)) stop("if provided, 'by' should be a (named) list of scalars")
  Z<-ogg$Z 
  nomeZ<-ogg$namesGZ$nameZ
  beta.noG<- fixef(ogg$lme.fit.noG) 
  beta.all<-fixef(ogg$lme.fit)
  #browser()
  
  beta.G<-beta.all[setdiff(names(fixef(ogg$lme.fit)), names(beta.noG))]
  
  
  #names(w)<-names(beta.G)<-sub("G\\.","",names(beta.G))
  #if(any(!names(by)%in%names(beta.G))) stop(" 'by'should include variable names affecting psi")
  #nomiCoef<-names(beta.noG)
  #browser()
  
  if(!is.null(by)) {
    #if(!all(sapply(by,is.numeric))){ #se c'e' qualche etichetta di categoriale
      a=by
      nomiOK<-const<-idList<-vector("list", length(a))
      value<-vector(,length(a))
      for(i in 1:length(a)) {
        nomeOK <- if(is.character(a[[i]])) paste(names(a[i]),a[[i]], sep="") else names(a[i])
        nomiOK[i]<- nomeOK <- paste("G", nomeOK, sep=".")
        bG<-beta.G[nomeOK]
        if(is.na(bG)) stop("the value supplied in 'by' is not coherent with the fit")
        if(!is.character(a[[i]])) bG <-bG*a[[i]]
        const[[i]]<-bG
        const[[i]]<- ifelse(is.na(const[[i]]),0,const[[i]])
        idList[[i]]<-names(bG)
        value[i]<-ifelse(is.character(a[[i]]),1,a[[i]])
        names(value[i])<-names(bG)
      }
      w<-c(1,value)
      final.names<-names(w)<-c("G0", unlist(nomiOK))
      beta.G<-beta.G[final.names]
      # const<-matrix(unlist(const),2, byrow=FALSE)
      # colnames(const)<-names(by)
      # nomiNOdiff <- names(which(colSums(const)==0))
      # if(length(nomiNOdiff)>0) warning("The value of", paste(" '", paste(nomiNOdiff, collapse=" and "),"' ",sep=""), 
      #                                    "supplied in 'by' does not modify the baseline slopes estimates", call. = FALSE)
      #   
      #   nomiCoef<- c(nomeZ, "U", unlist(idList))

    #} else {
      
      # w=c(1,rep(0, length(beta.G)-1))
      # names(w)<-final.names <-names(beta.G)
      # 
      # value<-unlist(by)
      # names(value)<- paste("G", names(value), sep=".")
      # w[names(value)]<-value
      
    #}
      #browser()
    r<-sum(beta.G*w)
    V<-V[final.names,final.names]
    SE.fit<-sqrt(drop(t(w)%*%V%*%w))
    zalpha<- -qnorm((1-conf.level)/2)
    ci<-cbind(r-zalpha*SE.fit, r+zalpha*SE.fit)
    #colnames(ci)<-paste("CI", "(", conf.level * 100,"%", ")", c(".low", ".up"), sep = "")
    colnames(ci)<-paste(" ", conf.level, c(".low", ".up"), sep = "")
    ris<-cbind("Est."=r, "St.Err"=SE.fit, ci)
    
  } else {
    r<-beta.G["G0"]
    SE.fit<-sqrt(V["G0","G0",drop=TRUE])
    zalpha<- -qnorm((1-conf.level)/2)
    ci<-cbind(r-zalpha*SE.fit, r+zalpha*SE.fit)
    #colnames(ci)<-paste("CI", "(", conf.level * 100,"%", ")", c(".low", ".up"), sep = "")
    colnames(ci)<-paste(" ", conf.level, c(".low", ".up"), sep = "")
    ris<-cbind("Est."=r, "St.Err"=SE.fit, ci)
  }
  if(use.link) {
    ris[,c(1,3,4)] <- ogg$misc$f.psi(ris[,c(1,3,4)])
    ris[,2]<-abs(ogg$misc$f.psi(ris[,1]))*ris[,2]
  }
  rownames(ris)<-"psi"
  
  #browser()
  
  ris
}

  
  




