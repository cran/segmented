`slope` <- function(ogg, parm, conf.level=0.95, rev.sgn=FALSE, APC=FALSE, .vcov=NULL, .coef=NULL, use.t=NULL,...,
      digits = max(4, getOption("digits") - 2)){

#   se e' un "newsegmented"
#    if(!is.null(ogg$R.slope)) {
#             covv<-old.coef.var(ogg)
#             ogg$coefficients<- covv$b
#             covv<-  covv$cov
#             ogg$psi<-old.psi(ogg)
#             ogg$nameUV<-old.nomi(ogg)
#             } else {
#             covv<-try(vcov(ogg,...), silent=TRUE)
#             }

        covv <- if(is.null(.vcov)) vcov(ogg, ...) else .vcov 
        estcoef<- if(is.null(.coef)) coef(ogg) else .coef
        
        if(length(estcoef)==0) stop("No coefficient in the object fit?")

        if(!all(dim(covv)==c(length(estcoef), length(estcoef)))) stop("dimension of cov matrix and estimated coeffs do not match", call. = FALSE)
        
        nomepsi<-rownames(ogg$psi) #OK
        nomeU<-ogg$nameUV$U
        nomeZ<-ogg$nameUV$Z
        if(missing(parm)) {
          nomeZ<- ogg$nameUV$Z
          if(length(rev.sgn)==1) rev.sgn<-rep(rev.sgn,length(nomeZ))
          } else {
            if(! all(parm %in% ogg$nameUV$Z)) {
              stop("invalid parm") } else {nomeZ<-parm}
            }
        if(length(rev.sgn)!=length(nomeZ)) rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        nomi<-names(estcoef)
        index<-vector(mode = "list", length = length(nomeZ))
        for(i in 1:length(nomeZ)) {
            index[[i]]<-match(c(nomeZ[i], names(ogg$indexU[[nomeZ[i]]])),  nomi)
        }
        if(!is.null(use.t)){
          k0<- if(use.t) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
        } else {
          k0<-if(inherits(ogg, "lm", which=TRUE)==2) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
        }
        #k<-if("lm"%in%class(ogg)) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
        Ris<-list()   
        #digits <- max(3, getOption("digits") - 3)
        rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        for(i in 1:length(index)){
            ind<-as.numeric(na.omit(unlist(index[[i]])))
            M<-matrix(1,length(ind),length(ind))
            M[row(M)<col(M)]<-0
            cof<-estcoef[ind]
            cof.out<-M%*%cof 

            if(!inherits(covv, "try-error")){
                cov.ok<-covv[ind,ind]
                cov.out<-M%*%cov.ok%*%t(M)
                se.out<-sqrt(diag(cov.out))
                k<-k0*se.out
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

