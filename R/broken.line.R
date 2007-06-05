`broken.line` <-
function(ogg,term=NULL,gap=FALSE,link=FALSE){
#returns the fitted straight lines from a `segmented' model
#term: a character vector meaning the segmented variable.
#     If NULL every segmented variable is considered and a matrix is returned
#gap: should the gap be accounted for? Currently unimplemented
#Problema: se ci sono poche osservazioni le linee risultano "brutte"..
        if(!"segmented"%in%class(ogg)) stop("A segmented model is needed")
        nomepsi<-rownames(ogg$psi) #OK
        nomeU<-ogg$nameUV[[1]]
        nomeZ<-ogg$nameUV[[3]]
        nomi<-names(coef(ogg))
        nomi<-nomi[-match(nomepsi,nomi)] #escludi i coef delle V
        index<-NULL
        for(i in 1:length(nomeZ)) index[[i]]<-grep(nomeZ[i], nomi, extended=FALSE)
        variabili<-Ris<-list()
        for(i in 1:length(index)){
            ind<-unlist(index[[i]])
            #M<-matrix(1,length(ind),length(ind))
            #M[row(M)<col(M)]<-0
            cof<-coef(ogg)[ind] #questo è il vettore di coef per la variabile seg
            Ris[[nomeZ[i]]]<-cof
            variabili[[nomeZ[i]]]<-data.matrix(ogg$mframe[names(cof)])
            }
        o<-mapply(function(xx,yy)drop(xx%*%yy),variabili,Ris)+coef(ogg)["(Intercept)"]
        if(!is.null(term)) o<-o[,term]
        if(inherits(ogg,what="glm",FALSE) && !link) o<-apply(o,2,ogg$family$linkinv)
        return(o)
        }
