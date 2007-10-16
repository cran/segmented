`ci.psi` <-
function(ogg,conf.level=0.95,Fieller=FALSE){
#restituisce CI per i psi
        if(!"segmented"%in%class(ogg)) stop("A segmented model is needed")
        nomeZ<-ogg$nameUV[[3]] #nomi delle variabili segmented
        Ris<-list()
        digits <- max(3, getOption("digits") - 3)
        for(i in 1:length(nomeZ)){
          id<-grep(nomeZ[i], rownames(ogg$psi), extended=FALSE)
          psi<-ogg$psi[id,2]
          se.psi<-ogg$psi[id,3] #SE of psi
          k<-abs(qnorm((1-conf.level)/2))*se.psi
          r<-cbind(psi,psi-k,psi+k)
          colnames(r)<-c("Est.",paste("CI","(",conf.level*100,"%",")",c(".l",".u"),sep=""))
          rownames(r)<-rownames(ogg$psi)[id]
          Ris[[nomeZ[i]]]<-signif(r,digits)
          }
        Ris
        }

