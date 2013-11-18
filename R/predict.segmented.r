predict.segmented<-function(object, newdata, ...){
#rev: 30/10/2013: it seems to work correctly, even with the minus variable (null right slope..)
#BUT problems if type="terms" (in realtà funziona, il problema è che
#     restituisce una colonna per "x", "U.x", "psi.x".. (Eventualmente si dovrebbero sommare..)
#returns predictions from a segmented fit
#object: the segmented fit
#newdata: it can be
#     1) an (unnamed) vector when there is a single segmented covariate;
#     2) a dataframe including *all* the variables in the model
#     3) a dataframe including a single segmented variable
#     4) missing, the model.frame is used.
  dummy.matrix<-function(x.values, x.name, obj.seg){
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
        V<-ifelse(newZ>PSI,-1,0)
        dummy2<- if(k==1) V*diffSlope  else V%*%diag(diffSlope) #t(diffSlope*t(-I(newZ>PSI)))

        newd<-cbind(x.values,dummy1,dummy2)
        colnames(newd)<-c(nameZ,nameU, nameV)
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
  nameU<-object$nameUV$U
  nameV<-object$nameUV$V
  nameZ<-object$nameUV$Z
  if(missing(newdata)){
    newd.ok<-model.frame(object)
  } else {
    if(is.vector(newdata)) {
      if(length(object$nameUV$Z)>1) stop("with multiple segmented variables, 'newdata' has to be a dataframe")
      newdata<-data.frame(newdata)
      names(newdata)<- object$nameUV$Z
      }
  if(!is.data.frame(newdata)) stop('newdata has to be a dataframe!')
  n<-nrow(newdata)
  newd.NOseg<-newdata[,setdiff(names(newdata), nameZ), drop=FALSE] #part of newdata including non-segmented variables only
  #what about factors????
  if(ncol(newd.NOseg)==ncol(newdata)) stop("'newdata' does not include any segmented variable")
  #se sei interessato soltanto ad una segmented variable
  if(ncol(newdata)<=1) nameZ<-intersect(nameZ, names(newdata))

  r<-NULL
  for(i in 1:length(nameZ)){
      x.name <- nameZ[i]
      x.values<-newdata[[x.name]]
      r[[i]]<-dummy.matrix(x.values, x.name, object)
      }
  newd.ok<-data.frame(matrix(unlist(r), nrow=n, byrow = FALSE))
  newd.ok<-cbind(newd.NOseg,newd.ok)
  names(newd.ok)<-sapply(r, colnames)
  #the dataframe to be passed to predict.lm() should include all the variables
  nomiRes<-setdiff(names(model.frame(object))[-1], names(newd.ok))

  NoSegLinMatrix<-matrix(0,nrow(newd.ok), length(nomiRes), dimnames=list(NULL, nomiRes))
  if(length(names(object$xlevels))>0){#se ci sono factor
          for(j in names(object$xlevels)){
              NoSegLinMatrix[,j]<-object$xlevels[[j]][1]
              NoSegLinMatrix[,j]<-factor(NoSegLinMatrix[,j])
              }
          }
  newd.ok <- cbind(newd.ok, NoSegLinMatrix)
  } #end of if(missing(newdata)) { } else {..}
  f<-if(inherits(object, what = "glm", which = FALSE)) predict.glm(object, newdata=newd.ok, ...) else predict.lm(object, newdata=newd.ok, ...)
  return(f)
  }
