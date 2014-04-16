broken.line<-function(ogg, term=NULL, link=TRUE, interc=TRUE, se.fit=TRUE){
#ogg: l'oggetto segmented
#term: una lista *nominata* con i valori rispetto a cui calcolare i fitted
#   OPPURE una stringa per indicare la variabile segmented OPPURE NULL (se c'è solo una variabile)
  dummy.matrix<-function(x.values, x.name, obj.seg, psi.est=TRUE){
    #given the segmented fit 'obj.seg' and a segmented variable x.name with corresponding values x.values,
    #this function simply returns a matrix with columns (x, (x-psi)_+, -b*I(x>psi))
    #or  ((x-psi)_+, -b*I(x>psi)) if obj.seg does not include the coef for the linear "x"
        n<-length(x.values)
        nameU<- grep(x.name, obj.seg$nameUV$U, value = TRUE)
        nameV<- grep(x.name, obj.seg$nameUV$V, value = TRUE)

        diffSlope<-coef(obj.seg)[nameU]
        est.psi<-obj.seg$psi[nameV,2]

        k<-length(est.psi)

        PSI <- matrix(rep(est.psi, rep(n, k)), ncol = k)
        newZ<-matrix(x.values, nrow=n,ncol=k, byrow = FALSE)

        dummy1<-pmax(newZ-PSI,0)
        
        if(psi.est){
          V<-ifelse(newZ>PSI,-1,0)
          dummy2<- if(k==1) V*diffSlope  else V%*%diag(diffSlope) #t(diffSlope*t(-I(newZ>PSI)))
          newd<-cbind(x.values,dummy1,dummy2)
          colnames(newd)<-c(x.name,nameU, nameV)
          } else {
          newd<-cbind(x.values,dummy1)
          colnames(newd)<-c(x.name,nameU)
          }
        if(!x.name%in%names(coef(obj.seg))) newd<-newd[,-1]
        return(newd)
  #  nullLeftSlope<-FALSE
  #  if(!nameZ%in%all.vars(object$orig.call)) nullLeftSlope<-TRUE
  #  if(ncol(newdata)>=2){
  #      newd.ok<-cbind(newd, newdata[,-match(nameZ,names(newdata)),drop=FALSE])
  #      names(newd.ok)<- c(nameZ,nameU, nameV,names(newdata)[-match(nameZ,names(newdata))])
  #      } else {
  #      newd.ok<-newd
  #      names(newd.ok)<-c(nameZ,nameU, nameV)
  #      }
  #  if(nullLeftSlope) newd.ok<-newd.ok[,-1]
    }
#--------------
    xvalues<-term
    nomeV <- ogg$nameUV$V
    nomeU <- ogg$nameUV$U
    nomeZ <- ogg$nameUV$Z
    n.seg<-length(nomeZ)
    if(is.null(xvalues)){
      if(n.seg>1) stop("there are multiple segmented covariates. Please specify which one.")
      xvalues<-ogg$model[nomeZ]
      }
    if(is.character(xvalues)){
          if(!xvalues %in% nomeZ) stop("'xvalues' is not a segmented covariate")
          xvalues<-ogg$model[xvalues]
      }
    nomeOK<-names(xvalues)
    if(!nomeOK %in% nomeZ) stop("'names(xvalues)' is not a segmented covariate")

    #if(n.seg>1 && !is.list(x.values)) stop("with multiple segmented covariates, please specify a named dataframe")
    #x.values<-data.frame(x.values)
    #names(x.values)<-nomeZ

    nomi <- names(coef(ogg))
    nomiSenzaV <- nomiSenzaU <- nomi
    nomiSenzaU[match(nomeU, nomi)] <- ""
    nomiSenzaV[match(nomeV, nomi)] <- ""

    index <- vector(mode = "list", length = length(nomeZ))
    for (i in 1:n.seg) {
        index[[i]] <- c(match(nomeZ[i], nomi),
            grep(paste("\\.", nomeZ[i], "$", sep = ""), nomiSenzaV, value = FALSE),
            grep(paste("\\.", nomeZ[i], "$", sep = ""), nomiSenzaU, value = FALSE))
            }
    ste.fit<-fit <- vector(mode = "list", length = length(nomeZ))
    for (i in 1:n.seg) {
        x.name <- nomeZ[i]
        X<-dummy.matrix(unlist(xvalues), x.name, ogg)
        ind <- as.numeric(na.omit(unlist(index[[i]])))
        if(interc && "(Intercept)"%in%nomi) {
          ind<- c(match("(Intercept)",nomi),ind)
          X<-cbind(1,X)
          }
        cof <- coef(ogg)[ind]
        fit[[i]]<-drop(X%*%cof)
        ste.fit[[i]] <- if(!se.fit) 10 else sqrt(rowSums((X %*% vcov(ogg)[ind,ind]) * X)) #sqrt(diag(X%*%Var%*%t(X)))
        }
        names(fit)<- names(ste.fit)<- nomeZ
        r<-list(fit=fit[[nomeOK]], se.fit=ste.fit[[nomeOK]])
        if (inherits(ogg, what = "glm", FALSE) && !link){
            r[[2]] <- ogg$family$mu.eta(r[[1]])*r[[2]]
            r[[1]] <- ogg$family$linkinv(r[[1]])
            }
        if(!se.fit) r<-r[1]
        return(r)
        }
