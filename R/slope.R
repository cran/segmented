`slope` <-
function(ogg, parm, conf.level=0.95, rev.sgn=FALSE, APC=FALSE, .vcov=NULL,...,
      digits = max(4, getOption("digits") - 2)){
#--

        f.U<-function(nomiU, term=NULL){
        #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
        #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".") #riscostruisce il nome con il "." (che era stato scomposto da strsplit())
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
#--        
#        if(!"segmented"%in%class(ogg)) stop("A 'segmented' model is requested")

        #commentato il 28/05/19. vcov() fa gia' questo controllo..
        # if(var.diff && length(ogg$nameUV$Z)>1) {
        #     var.diff<-FALSE
        #     warning("var.diff set to FALSE with multiple segmented variables", call.=FALSE)
        #     }

            #se e' un "newsegmented"
#    if(!is.null(ogg$R.slope)) {
#             covv<-old.coef.var(ogg)
#             ogg$coefficients<- covv$b
#             covv<-  covv$cov
#             ogg$psi<-old.psi(ogg)
#             ogg$nameUV<-old.nomi(ogg)
#             } else {
             covv<-try(vcov(ogg,...), silent=TRUE)
#             }

        covv <- if(is.null(.vcov)) vcov(ogg, ...) else .vcov 
        if(!all(dim(covv)==c(length(coef(ogg)), length(coef(ogg))))) stop("Incorrect dimension of cov matrix", call. = FALSE)
             
        nomepsi<-rownames(ogg$psi) #OK
        nomeU<-ogg$nameUV$U
        nomeZ<-ogg$nameUV$Z
        if(missing(parm)) {
          nomeZ<- ogg$nameUV$Z
          if(length(rev.sgn)==1) rev.sgn<-rep(rev.sgn,length(nomeZ))
          }
             else {
                if(! all(parm %in% ogg$nameUV$Z)) {stop("invalid parm")}
                  else {nomeZ<-parm}
                  }
        if(length(rev.sgn)!=length(nomeZ)) rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        nomi<-names(coef(ogg))
        nomi<-nomi[-match(nomepsi,nomi)] #escludi i coef delle V
        index<-vector(mode = "list", length = length(nomeZ))
        for(i in 1:length(nomeZ)) {
#--->             DA RIMUOVERE E SOSTITUIRE CON QUELLI DI SUBITO DOPO?
#            #id.cof.U<-grep(paste("\\.",nomeZ[i],"$",sep=""), nomi, value=FALSE)
#            #psii<-ogg$psi[grep(paste("\\.",nomeZ[i],"$",sep=""), rownames(ogg$psi), value=FALSE),2]
#            #id.cof.U<- match(grep(nomeZ[i],   ogg$nameUV$U, value=TRUE), nomi)
#            #psii<-ogg$psi[grep(nomeZ[i],   ogg$nameUV$V, value=TRUE),2]
#            #il paste con "$" (paste("\\.",nomeZ[i],"$",sep="")) e' utile perche' serve a distinguere variabili con nomi simili (ad es., "x" e "xx")
#            #Comunque nella versione dopo la 0.3-1.0 ho (FINALMENTE) risolto mettendo f.U
#            id.cof.U<- f.U(ogg$nameUV$U, nomeZ[i])
#            #id.cof.U e' la posizione nel vettore ogg$nameUV$U; la seguente corregge per eventuali variabili che ci sono prima (ad es., interc)
#            id.cof.U<- id.cof.U + (match(ogg$nameUV$U[1], nomi)-1)            
#            psii<- ogg$psi[f.U(ogg$nameUV$V, nomeZ[i]) , "Est."]
#            id.cof.U <- id.cof.U[order(psii)]            
#--->            
            
            #questi funzionano anche con oggetti da segreg
            nomiPsi<-grep(paste(".", nomeZ[i], sep="") , ogg$nameUV$V, value=TRUE)
            psii<- ogg$psi[nomiPsi , "Est."] 
            nomiU<-grep(paste(".", nomeZ[i], sep="") , ogg$nameUV$U, value=TRUE)
            #cof<-coef(ogg)[nomiU]
            id.cof.U<- match(nomiU, names(coef(ogg))) #prima era names(ogg$coefficients)
            index[[i]]<-c(match(nomeZ[i],nomi), id.cof.U)
            }
        Ris<-list()   
        #digits <- max(3, getOption("digits") - 3)
        rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        
#         transf=c("x","1")
#        if( (length(transf)!=2) || !(length(transf)==1 && transf=="APC")) stop("'error in transf'")
#        if(transf=="APC") transf<-c("100*(exp(x)-1)", "100*exp(x)")
#        my.f<-function(x)eval(parse(text=transf[1]))
#        my.f.deriv<-function(x)eval(parse(text=transf[2]))
       

        #browser()
                  
        for(i in 1:length(index)){
            ind<-as.numeric(na.omit(unlist(index[[i]])))
            M<-matrix(1,length(ind),length(ind))
            M[row(M)<col(M)]<-0
            cof<-coef(ogg)[ind]
            cof.out<-M%*%cof 

            if(!inherits(covv, "try-error")){
            #if(class(covv)[1]!="try-error"){ #inherits(covv, "try-error") #is(covv, "try-error")
                cov.ok<-covv[ind,ind]
                cov.out<-M%*%cov.ok%*%t(M)
                se.out<-sqrt(diag(cov.out))
                k<-if("lm"%in%class(ogg)) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
                k<-k*se.out
                ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
                } else {
                ris<-cbind(cof.out, NA, NA, NA,NA)                
                }
                cin<-paste("CI","(",conf.level*100,"%",")",c(".l",".u"),sep="")
            #se la left slope e' nulla....
            #if(identical(length(ind),length(grep(paste("\\.",nomeZ[i],"$",sep=""), nomeU)))){
            if(!nomeZ[i]%in%nomi){            
                    ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)
                    }
            if(rev.sgn[i]){
                ris<-cbind(-ris[nrow(ris):1,1],ris[nrow(ris):1,2],-ris[nrow(ris):1,3],
                      -ris[nrow(ris):1,5],-ris[nrow(ris):1,4])
                      }
            nomeT<-if("lm"%in%class(ogg)) "t value" else "z value"
            dimnames(ris)<-list(paste("slope", 1:nrow(ris), sep=""),c("Est.","St.Err.",nomeT,cin[1],cin[2]))
            if(APC) ris<-100*(exp(ris[,c(1,4,5)])-1)
            Ris[[nomeZ[i]]]<- signif(ris,digits)
            
                } #end loop i
            #if(!missing(parm)){
            #    if(!all(parm %in% ogg$nameUV[[3]])) stop("invalid parm") else Ris<-Ris[parm]}
            Ris
            }

