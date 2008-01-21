`slope` <-
function(ogg, parm, conf.level=0.95, rev.sgn=FALSE){
#Ad ogni chiamata di grep() è stato aggiunto l'argomento extended=F per consentire Z=log(x)
        if(!"segmented"%in%class(ogg)) stop("A segmented model is needed")
        nomepsi<-rownames(ogg$psi) #OK
        nomeU<-ogg$nameUV[[1]]
        nomeZ<-ogg$nameUV[[3]]
        if(missing(parm)) {
          nomeZ<- ogg$nameUV[[3]]
          if(length(rev.sgn)==1) rev.sgn<-rep(rev.sgn,length(nomeZ))
          }
             else {
                if(! all(parm %in% ogg$nameUV[[3]])) {stop("invalid parm")}
                  else {nomeZ<-parm}
                  }
        if(length(rev.sgn)!=length(nomeZ)) rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        nomi<-names(coef(ogg))
        nomi<-nomi[-match(nomepsi,nomi)] #escludi i coef delle V
        index<-vector(mode = "list", length = length(nomeZ))
        for(i in 1:length(nomeZ)) index[[i]]<-grep(nomeZ[i], nomi, extended=FALSE) 
        Ris<-list()   
        digits <- max(3, getOption("digits") - 3)
        rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        for(i in 1:length(index)){
            ind<-unlist(index[[i]])
            M<-matrix(1,length(ind),length(ind))
            M[row(M)<col(M)]<-0
            cof<-coef(ogg)[ind]
            covv<-vcov(ogg)[ind,ind] 
            cof.out<-M%*%cof 
            cov.out<-M%*%covv%*%t(M)
            se.out<-sqrt(diag(cov.out))
            k<-abs(qnorm((1-conf.level)/2))*se.out
            ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
            cin<-paste("CI","(",conf.level*100,"%",")",c(".l",".u"),sep="")
            #se la left slope è nulla....
            if(identical(length(ind),length(grep(nomeZ[i], nomeU, extended=FALSE)))){
                    ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)}
            if(rev.sgn[i]){
                ris<-cbind(-ris[nrow(ris):1,1],ris[nrow(ris):1,2],-ris[nrow(ris):1,3],
                      -ris[nrow(ris):1,5],-ris[nrow(ris):1,4])}
            dimnames(ris)<-list(paste("slope", 1:nrow(ris), sep=""),c("Est.","St.Err.","t value",cin[1],cin[2]))
            Ris[[nomeZ[i]]]<-signif(ris,digits)
                } #end loop i
            #if(!missing(parm)){
            #    if(!all(parm %in% ogg$nameUV[[3]])) stop("invalid parm") else Ris<-Ris[parm]}
            Ris
            }

