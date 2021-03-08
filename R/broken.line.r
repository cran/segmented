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
        estcoef <- if(is.null(.coef)) coef(obj.seg) else .coef
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
        if(length(isV)==1) isV<-c(FALSE,isV)
        n<-length(x.values)
        #le seguenti righe selezionavano (ERRONEAMENTE) sia "U1.x" sia "U1.neg.x" (se "x" e "neg.x" erano segmented covariates)
        #nameU<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$U, value = TRUE)
        #nameV<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$V, value = TRUE)
        nameU<-obj.seg$nameUV$U[f.U(obj.seg$nameUV$U,x.name)]
        nameV<-obj.seg$nameUV$V[f.U(obj.seg$nameUV$V,x.name)]

        diffSlope<-estcoef[nameU]
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
          newd<-cbind(x.values,dummy1,dummy2)
          colnames(newd)<-c(x.name,nameU, nameV)
          } else {
          newd<-cbind(x.values,dummy1)
          colnames(newd)<-c(x.name,nameU)
          }
        if(!x.name%in%names(estcoef)) newd<-newd[,-1,drop=FALSE]
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
        return(newd)
    } #end dummy.matrix()
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
    xvalues<-term
    nomeV <- ogg$nameUV$V
    nomeU <- ogg$nameUV$U
    nomeZ <- ogg$nameUV$Z
    n.seg<-length(nomeZ)
    if(is.null(xvalues)){
      if(n.seg>1) stop("there are multiple segmented covariates. Please specify one.")
      xvalues<-ogg$model[nomeZ]
      }
    if(is.character(xvalues)){
          if(!xvalues %in% nomeZ) stop("'xvalues' is not a segmented covariate")
          xvalues<-ogg$model[xvalues]
      }
    nomeOK<-names(xvalues)
    if(length(nomeOK)>1) stop("Please specify one variable")
    if(!nomeOK %in% nomeZ) stop("'names(xvalues)' is not a segmented covariate")

    nomi <- names(estcoef)
    nomiSenzaV <- nomiSenzaU <- nomi
    nomiSenzaU[match(nomeU, nomi)] <- ""
    nomiSenzaV[match(nomeV, nomi)] <- ""
    index <- vector(mode = "list", length = length(nomeZ))
    for (i in 1:n.seg) {
        index[[i]] <- c(match(nomeZ[i], nomi),
            f.U(ogg$nameUV$U, nomeZ[i]) + (match(ogg$nameUV$U[1], nomi)-1),
            f.U(ogg$nameUV$V, nomeZ[i]) + (match(ogg$nameUV$V[1], nomi)-1))
            #match(setdiff(names(ogg$indexU[[x.name]]), )))
            }
    ste.fit<-fit <- vector(mode = "list", length = length(nomeZ))
    for(i in 1:n.seg) {
        x.name <- nomeZ[i]
        Xfit<-dummy.matrix(unlist(xvalues), x.name, ogg, isV=FALSE, .coef=estcoef)
        if(se.fit) X<-dummy.matrix(unlist(xvalues), x.name, ogg, isV=isV, .coef=estcoef)#<--NB: xvalues non varia con i!!! perche' farlo calcolare comunque? 
        ind <- as.numeric(na.omit(unlist(index[[i]])))
        
        idInterc<-grep("ntercept",names(estcoef))
        #if(idInterc!= grep("ntercept",rownames(covv))) stop("intercept name in coeff and vcov do not match")
         
        
        
        #if(interc && "(Intercept)"%in%nomi) {
        if(interc && length(idInterc)==1) {
          ind<- c(idInterc, ind)
          Xfit<-cbind(1,Xfit)
          if(se.fit) X<-cbind(1,X)
          }
        ind<-union(ind, match(names(ogg$indexU[[x.name]]), nomi))
#browser()
        cof <- estcoef[ind]
        fit[[i]]<-drop(Xfit%*%cof)
        if(se.fit) ste.fit[[i]] <- sqrt(rowSums((X %*% covv[ind,ind]) * X)) #sqrt(diag(X%*%Var%*%t(X)))
        #ste.fit[[i]] <- if(!se.fit) 10 else sqrt(rowSums((X %*% vcov.segmented(ogg,...)[ind,ind]) * X)) #sqrt(diag(X%*%Var%*%t(X)))
        }
        names(fit)<- names(ste.fit)<- nomeZ
        r<-list(fit=fit[[nomeOK]], se.fit=ste.fit[[nomeOK]])
        if (inherits(ogg, what = "glm", FALSE) && !link){
          if(se.fit) r[[2]] <- ogg$family$mu.eta(r[[1]])*r[[2]]
          r[[1]] <- ogg$family$linkinv(r[[1]])
            }
        if(!se.fit) r<-r[1] #metti r[[1]] se vuoi fare restituire un vettore
        return(r)
        }
