`intercept` <-
function(ogg, parm, rev.sgn=FALSE){#conf.level=0.95, var.diff=FALSE
        if(!"segmented"%in%class(ogg)) stop("A segmented model is needed")
        if(var.diff && length(ogg$nameUV$Z)>1) {
            var.diff<-FALSE
            warning("var.diff set to FALSE with multiple segmented variables", call.=FALSE)
            }
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
        Allpsi<-index<-vector(mode = "list", length = length(nomeZ))

        for(i in 1:length(nomeZ)) {
            id.cof.U<-grep(paste("\\.",nomeZ[i],"$",sep=""), nomi, value=FALSE)
            psii<-ogg$psi[grep(paste("\\.",nomeZ[i],"$",sep=""), rownames(ogg$psi), value=FALSE),2]
            Allpsi[[i]]<-psii
            id.cof.U <- id.cof.U[order(psii)]
            index[[i]]<-id.cof.U #solo diffSlope
            #index[[i]]<-c(match(nomeZ[i],nomi), id.cof.U) #id della left slope+ diffSlope
            }
        Ris<-list()
        digits <- max(3, getOption("digits") - 3)
        rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        alpha0<-alpha00<-coef(ogg)["(Intercept)"]


        for(i in 1:length(index)){
            ind<-as.numeric(na.omit(unlist(index[[i]])))
#            M<-matrix(1,length(ind),length(ind))
#            M[row(M)<col(M)]<-0
            cof<-coef(ogg)[ind]
            alpha<-vector(length=length(ind))
            for(j in 1:length(cof)){
                alpha[j]<-alpha0-Allpsi[[i]][j]*cof[j]
                alpha0<-alpha[j]
                }
            cof.out<-c(alpha00,alpha)
            ris<-matrix(cof.out)
            dimnames(ris)<-list(paste("intercept", 1:nrow(ris), sep=""),"Est.") #,"St.Err.",nomeT,cin[1],cin[2]))
#            covv<-vcov(ogg,var.diff=var.diff)[ind,ind]
#            cof.out<-M%*%cof
#            cov.out<-M%*%covv%*%t(M)
#            se.out<-sqrt(diag(cov.out))
#            k<-if("lm"%in%class(ogg)) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
#            k<-k*se.out
#            ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
#            cin<-paste("CI","(",conf.level*100,"%",")",c(".l",".u"),sep="")
#            #se la left slope è nulla....
#            if(identical(length(ind),length(grep(paste("\\.",nomeZ[i],"$",sep=""), nomeU)))){
#                    ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)}
#            if(rev.sgn[i]){
#                ris<-cbind(-ris[nrow(ris):1,1],ris[nrow(ris):1,2],-ris[nrow(ris):1,3],
#                      -ris[nrow(ris):1,5],-ris[nrow(ris):1,4])}
#            nomeT<-if("lm"%in%class(ogg)) "t value" else "z value"
#            dimnames(ris)<-list(paste("slope", 1:nrow(ris), sep=""),c("Est.","St.Err.",nomeT,cin[1],cin[2]))
            Ris[[nomeZ[i]]]<-signif(ris,digits)
                } #end loop i
            Ris
            }

