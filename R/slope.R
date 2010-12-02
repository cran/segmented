`slope` <-
function(ogg, parm, conf.level=0.95, rev.sgn=FALSE, var.diff=FALSE, APC=FALSE){
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
        index<-vector(mode = "list", length = length(nomeZ))
        for(i in 1:length(nomeZ)) {
            id.cof.U<-grep(paste("\\.",nomeZ[i],"$",sep=""), nomi, value=FALSE)
            psii<-ogg$psi[grep(paste("\\.",nomeZ[i],"$",sep=""), rownames(ogg$psi), value=FALSE),2]
            id.cof.U <- id.cof.U[order(psii)]            
            index[[i]]<-c(match(nomeZ[i],nomi), id.cof.U)
            }
        Ris<-list()   
        digits <- max(3, getOption("digits") - 3)
        rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        
#         transf=c("x","1")
#        if( (length(transf)!=2) || !(length(transf)==1 && transf=="APC")) stop("'error in transf'")
#        if(transf=="APC") transf<-c("100*(exp(x)-1)", "100*exp(x)")
#        my.f<-function(x)eval(parse(text=transf[1]))
#        my.f.deriv<-function(x)eval(parse(text=transf[2]))
          
        for(i in 1:length(index)){
            ind<-as.numeric(na.omit(unlist(index[[i]])))
            M<-matrix(1,length(ind),length(ind))
            M[row(M)<col(M)]<-0
            cof<-coef(ogg)[ind]
            covv<-vcov(ogg,var.diff=var.diff)[ind,ind]
            cof.out<-M%*%cof 
            cov.out<-M%*%covv%*%t(M)
            se.out<-sqrt(diag(cov.out))
            k<-if("lm"%in%class(ogg)) abs(qt((1-conf.level)/2,df=ogg$df.residual)) else abs(qnorm((1-conf.level)/2))
            k<-k*se.out
            ris<-cbind(cof.out,se.out,(cof.out/se.out),(cof.out-k),(cof.out+k))
            cin<-paste("CI","(",conf.level*100,"%",")",c(".l",".u"),sep="")
            #se la left slope è nulla....
            if(identical(length(ind),length(grep(paste("\\.",nomeZ[i],"$",sep=""), nomeU)))){
                    ris<-rbind(c(0,rep(NA,(ncol(ris)-1))),ris)}
            if(rev.sgn[i]){
                ris<-cbind(-ris[nrow(ris):1,1],ris[nrow(ris):1,2],-ris[nrow(ris):1,3],
                      -ris[nrow(ris):1,5],-ris[nrow(ris):1,4])}
            nomeT<-if("lm"%in%class(ogg)) "t value" else "z value"
            dimnames(ris)<-list(paste("slope", 1:nrow(ris), sep=""),c("Est.","St.Err.",nomeT,cin[1],cin[2]))
            if(APC) ris<-100*(exp(ris[,c(1,4,5)])-1)
            Ris[[nomeZ[i]]]<-signif(ris,digits)
            
                } #end loop i
            #if(!missing(parm)){
            #    if(!all(parm %in% ogg$nameUV[[3]])) stop("invalid parm") else Ris<-Ris[parm]}
            Ris
            }

