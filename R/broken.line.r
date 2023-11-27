broken.line<-function(ogg, term=NULL, link=TRUE, interc=TRUE, se.fit=TRUE, isV=FALSE, .vcov=NULL, .coef=NULL, ...){
#ogg: l'oggetto segmented
#term: una lista *nominata* con i valori rispetto a cui calcolare i fitted
#   OPPURE una stringa per indicare la variabile segmented OPPURE NULL (se c'e' solo una variabile)
#is: 2 valori T/F per indicare se le variabili U e V nella matrice X, andrebbero sostituite con le versioni ind-smooth prima di calcolare var(X\hat\beta)
#...: argomenti da passare a vcov.segmented(): per esempio var.diff, is, p.df  
    if(length(isV)==1) isV<-c(FALSE,isV)
    dummy.matrix<-NULL
    dummy.matrix<-function(x.values, x.name, obj.seg, psi.est=TRUE, isV=FALSE, .coef=NULL){
        #given the segmented fit 'obj.seg' and a segmented variable x.name with corresponding values x.values,
        #this function simply returns a matrix with columns (x, (x-psi)_+, -b*I(x>psi))
        #or  ((x-psi)_+, -b*I(x>psi)) if obj.seg does not include the coef for the linear "x"
        f.U<-function(nomiU, term=NULL){
          #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
         #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
        estcoef <- if(is.null(.coef)) coef(obj.seg) else .coef
        if(length(isV)==1) isV<-c(FALSE,isV)
        n<-length(x.values)
        #le seguenti righe selezionavano (ERRONEAMENTE) sia "U1.x" sia "U1.neg.x" (se "x" e "neg.x" erano segmented covariates)
        #nameU<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$U, value = TRUE)
        #nameV<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$V, value = TRUE)
        nameU<-obj.seg$nameUV$U[f.U(obj.seg$nameUV$U,x.name)]
        nameV<-obj.seg$nameUV$V[f.U(obj.seg$nameUV$V,x.name)]

        #browser()
        if(is.null(obj.seg$constr)){
          diffSlope<-estcoef[nameU]
        } else {
          diffSlope<-drop(obj.seg$constr$invA.RList[[match(x.name, obj.seg$nameUV$Z)]]%*%estcoef[nameU])[-1]
        }
        
        est.psi<-obj.seg$psi[nameV, "Est."]
        se.psi<-obj.seg$psi[nameV, "St.Err"]
        if(any(is.na(se.psi))) stop("The St.Err. of psi is NA", call. = FALSE)
        k<-length(est.psi)

        PSI <- matrix(rep(est.psi, rep(n, k)), ncol = k)
        SE.PSI <- matrix(rep(se.psi, rep(n, k)), ncol = k)
        newZ<-matrix(x.values, nrow=n,ncol=k, byrow = FALSE)

        dummy1<-if(isV[1]) (newZ-PSI)*pnorm((newZ-PSI)/SE.PSI) else (newZ-PSI)*(newZ>PSI) #pmax(newZ-PSI,0)
        
        if(psi.est){
          V<-if(isV[2]) -pnorm((newZ-PSI)/SE.PSI) else -(newZ>PSI) #ifelse(newZ>PSI,-1,0)
          dummy2<- if(k==1) V*diffSlope  else V%*%diag(diffSlope) #t(diffSlope*t(-I(newZ>PSI)))
          colnames(dummy2)<- nameV
          newd<-cbind(x.values,dummy1,dummy2)
          colnames(newd)<-c(x.name,sub("psi","U", nameV), nameV)
          #colnames(newd)[1]<- x.name
          #colnames(newd)<-c(x.name,nameU, nameV)
          } else {
          newd<-cbind(x.values,dummy1)
          colnames(newd)<-c(x.name, sub("psi","U", nameV))
          #colnames(newd)[1]<- x.name
          #colnames(newd)<-c(x.name,nameU)
          }
        #if(!x.name%in%names(estcoef)) newd<-newd[,-1,drop=FALSE]
        #aggiungi (eventualmente) le colonne relative ai psi noti
        all.psi<-obj.seg$indexU[[x.name]]
        if(length(all.psi)!=k){
          nomi.psi.noti<-setdiff(names(all.psi),nameU)
          psi.noti<-setdiff(all.psi, est.psi)
          PSI.noti <- matrix(rep(psi.noti, rep(n, length(psi.noti))), ncol = length(psi.noti))
          nomi<-c(colnames(newd),nomi.psi.noti)
          newZ<-matrix(newZ, nrow=nrow(newZ), ncol=length(psi.noti))
          newd<-cbind(newd, (newZ-PSI.noti)*(newZ>PSI.noti))
          colnames(newd)<-nomi
        }
        rownames(newd)<-NULL
        return(newd)
    } #end dummy.matrix()
    #--------------
    blockdiag <- function(...) {
      args <- list(...)
      nc <- sapply(args,ncol)
      cumnc <- cumsum(nc)
      ##  nr <- sapply(args,nrow)
      ## NR <- sum(nr)
      NC <- sum(nc)
      rowfun <- function(m,zbefore,zafter) {
        cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
              matrix(0,ncol=zafter,nrow=nrow(m)))
      }
      ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
      for (i in 2:length(args)) {
        ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
      }
      ret
    }
#--------------
    f.U<-function(nomiU, term=NULL){
        #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
        #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
#-------------
    estcoef<- if(is.null(.coef)) coef(ogg) else .coef
    if(se.fit) {
      covv<- if(is.null(.vcov)) vcov.segmented(ogg, ...) else .vcov
      #---Dalla versione 1.2.0 (20/06/20) ho eliminato il controllo sotto per consentire l'utilizzo
      #----di modelli che restituivano una cov con dimensione diversa dal numero dei coeff lineari (ad es., censReg)
      #if(!all(dim(.vcov)==c(length(ogg$coef), length(ogg$coef)))) stop("Incorrect dimension of cov matrix", call. = FALSE)
      if(!all(dim(covv)==c(length(estcoef), length(estcoef)))) stop("dimension of cov matrix and estimated coeffs do not match", call. = FALSE)
    }
    #browser()
    nomeZ <- ogg$nameUV$Z
    if(is.null(term)){
      term <- nomeZ[1]
      xvalues<-ogg$model[term]
    } else {
      if(is.character(term)) term<- ogg$model[term]
      if(!is.list(term)) stop("term should be a named list")
      if(!names(term)%in%nomeZ) stop("term is not a segmented variable")
      xvalues<-term
      term<-names(term)
    }

    n.seg<-1
    # if(is.null(xvalues)){
    #   if(n.seg>1) stop("there are multiple segmented covariates. Please specify one.")
    #   xvalues<-ogg$model[nomeZ]
    #   }
    # if(is.character(xvalues)){
    #       if(!xvalues %in% nomeZ) stop("'xvalues' is not a segmented covariate")
    #       xvalues<-ogg$model[xvalues]
    #   }
    # nomeOK<-names(xvalues)
    # if(length(nomeOK)>1) stop("Please specify one variable")
    # if(!nomeOK %in% nomeZ) stop("'names(xvalues)' is not a segmented covariate")

    nomi <- names(estcoef)
    #nomiSenzaV <- nomiSenzaU <- nomi
    #nomiSenzaU[match(nomeU, nomi)] <- ""
    #nomiSenzaV[match(nomeV, nomi)] <- ""
    
    idInterc<-grep("ntercept",names(estcoef))
    
    ste.fit<-fit <- vector(mode = "list", length = 1)
    n.seg<-1
    #browser()
    
    ind <- match(c(term,ogg$nameUV$U[grep(term, ogg$nameUV$U)]), nomi, 0)
    ind<-ind[ind!=0]
    indV <- match(c(ogg$nameUV$V[grep(term, ogg$nameUV$V)]), nomi, 0)
    #Xfit<-dummy.matrix(unlist(xvalues), term, ogg, isV=FALSE, .coef=estcoef)
    #if(se.fit) 
    X<-dummy.matrix(unlist(xvalues), term, ogg, isV=isV, .coef=estcoef)
     #browser() 
    if(is.null(ogg$constr)){
        cof <- estcoef[ind]
        if(!term%in%names(estcoef)) X <-X[,-1,drop=FALSE]
        } else {
          #idU.i <- match(nomeU.i, names(estcoef))
          cof<-drop(ogg$constr$invA.RList[[match(term, ogg$nameUV$Z)]]%*%estcoef[ind])
          names(cof)<-c(term, paste("U",1:(length(cof)-1),".",term,sep="" ))
          #estcoef<-append(estcoef[-ind], cof, after=ind[1]-1)
        }
    
    idV <- match(grep(term, ogg$nameUV$V, value = TRUE),colnames(X))
    Xfit <- X[, -idV, drop=FALSE]
      
    fit<-drop(Xfit%*%cof)
    if(interc && length(idInterc)==1){
      fit <- fit + estcoef[idInterc]
      #if(se.fit) 
      X<-cbind(1,X)
      ind <- c(idInterc, ind)
    }

    #browser()
    
    if(se.fit) {
      V <- covv[c(ind, indV), c(ind, indV)]
      if(!is.null(ogg$constr)){
        B <- ogg$constr$invA.RList[[match(term, ogg$nameUV$Z)]]
        B <- do.call(blockdiag, list(diag(interc), B, diag(length(indV)))) 
        V <- B %*% V %*% t(B)
      } 
      #else {
      #  V <- vcov(object) 
      #  X <- X[,colnames(V)]
      #}
      
      ste.fit <- sqrt(rowSums((X %*% V) * X)) #sqrt(diag(X%*%Var%*%t(X)))
      if(inherits(ogg, what = "glm", FALSE) && !link) {
        ste.fit <- ogg$family$mu.eta(fit)*ste.fit
        fit <- ogg$family$linkinv(fit)
      }
      fit<-list(fit=fit, se.fit=ste.fit)
      } else {
        if(inherits(ogg, what = "glm", FALSE) && !link) fit <- ogg$family$linkinv(fit)
        fit<-list(fit=fit)
      }
    return(fit)
}


