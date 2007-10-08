`slope` <-
function(ogg, level=0.95){
#2/10/03; 21/10/03; 4/11/03; 24/02/04; 25/02/04; 28/04/06
#Ad ogni chiamata di grep() è stato aggiunto l'argomento extended=F per consentire Z=log(x)
#returns Est., St.Err., t value and Conf Interv (1-level) for the slopes
#      of the variables having a segmented relationship in the model.
#ogg: a segmented object returned by segmented().
        if(!"segmented"%in%class(ogg)) stop("A segmented model is needed")
        nomepsi<-rownames(ogg$psi) #OK
        nomeU<-ogg$nameUV[[1]]
        nomeZ<-ogg$nameUV[[3]]
        nomi<-names(coef(ogg))
        nomi<-nomi[-match(nomepsi,nomi)] #escludi i coef delle V
        index<-vector(mode = "list", length = length(nomeZ))
        for(i in 1:length(nomeZ)) index[[i]]<-grep(nomeZ[i], nomi, extended=FALSE) 
        Ris<-list()   
        digits <- max(3, getOption("digits") - 3)
        for(i in 1:length(index)){
            ind<-unlist(index[[i]])
            M<-matrix(1,length(ind),length(ind))
            M[row(M)<col(M)]<-0
            cof<-coef(ogg)[ind]
            covv<-vcov(ogg)[ind,ind] 
            cof.out<-M%*%cof 
            cov.out<-M%*%covv%*%t(M)
            se.out<-sqrt(diag(cov.out))
            k<-abs(qnorm((1-level)/2))*se.out
            ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
            cin<-paste("CI","(",level*100,"%",")",c(".l",".u"),sep="")
            #se la left slope è nulla....
            if(identical(length(ind),length(grep(nomeZ[i], nomeU, extended=FALSE)))){
                    ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)}
            dimnames(ris)<-list(paste("slope", 1:nrow(ris), sep=""),c("Est.","St.Err.","t value",cin[1],cin[2]))
            Ris[[nomeZ[i]]]<-signif(ris,digits)
                }
            Ris
            }

