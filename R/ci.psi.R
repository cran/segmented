`ci.psi` <-
function(ogg,level=0.95,Fieller=FALSE){
#restituisce CI per i psi
        if(!"segmented"%in%class(ogg)) stop("A segmented model is needed")
        nomeZ<-ogg$nameUV[[3]] #nomi delle variabili segmented
        Ris<-list()
        for(i in 1:length(nomeZ)){
          id<-grep(nomeZ[i], rownames(ogg$psi), extended=FALSE)
          psi<-ogg$psi[id,2]
          se.psi<-ogg$psi[id,3] #SE of psi
          k<-abs(qnorm((1-level)/2))*se.psi
          r<-cbind(psi,psi-k,psi+k)
          colnames(r)<-c("Est.",paste("CI","(",level*100,"%",")",c(".l",".u"),sep=""))
          rownames(r)<-rownames(ogg$psi)[id]
          Ris[[nomeZ[i]]]<-r
          }
        Ris
        }

