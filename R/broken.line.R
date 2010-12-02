`broken.line` <-
function(ogg,term=NULL,gap=FALSE,linkinv=FALSE,interc=TRUE){
#returns the fitted values according to each piecewise relationship from a `segmented' model
#term: a character vector meaning the segmented variable.
#     If NULL every segmented variable is considered and a matrix is returned
#gap: should the gap be accounted for? 
#Problema: se ci sono poche osservazioni le linee risultano "brutte"..
#1/12/10
        if(!"segmented"%in%class(ogg)) stop("A segmented model is requested")
        nomepsi<-rownames(ogg$psi) #OK
        nomeU<-ogg$nameUV$U
        nomeZ<-ogg$nameUV$Z
        nomiSenzaV<-nomiSenzaU<-nomi<-names(coef(ogg))
#        nomiSenzaV<-nomi[-match(nomepsi,nomi)] #setdiff(nomi, nomepsi)
        nomiSenzaU[match(nomeU,nomi)]<-""
        nomiSenzaV[match(nomepsi,nomi)]<-""
        index<-vector(mode = "list", length = length(nomeZ))
        for(i in 1:length(nomeZ)) {
          index[[i]]<-c( match(nomeZ[i],nomi), 
                    grep(paste("\\.",nomeZ[i],"$",sep=""), nomiSenzaV,value=FALSE))
          if(gap) index[[i]]<-c(index[[i]], 
                      grep(paste("\\.",nomeZ[i],"$",sep=""), nomiSenzaU,value=FALSE)
                    )
          }
        variabili<-Ris<-list()
        for(i in 1:length(index)){
            ind<-as.numeric(na.omit(unlist(index[[i]])))
            #M<-matrix(1,length(ind),length(ind))
            #M[row(M)<col(M)]<-0
            cof<-coef(ogg)[ind] #questo è il vettore di coef per la variabile seg
            Ris[[nomeZ[i]]]<-cof
            variabili[[nomeZ[i]]]<-data.matrix(ogg$model[names(cof)])
            }
        ris<-mapply(function(xx,yy)drop(xx%*%yy),variabili,Ris)
        if(interc) ris<-ris + coef(ogg)["(Intercept)"]
        if(!is.null(term)) ris<-ris[,term]
        if(inherits(ogg,what="glm",FALSE) && linkinv) ris<-apply(ris,2,ogg$family$linkinv)
        return(ris)
        }


    