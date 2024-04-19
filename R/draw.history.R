draw.history<-function(obj,term,...){
#show.history() se c'e' stato boot restart potrebbe produrre un grafico 2x1 di "dev vs it" and "no.of distinct vs it"
#--
        f.U<-function(nomiU, term=NULL){
        #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
        #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
#--        
      if(missing(term)){
          if(length(obj$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-obj$nameUV$Z}
      }
       #browser() 
      opz<-list(...)
      range.ok<-obj$rangeZ[,term]
      id.ok<- f.U(rownames(obj$psi), term)
      est.psi<-obj$psi[id.ok,"Est."]
      if(is.null(opz$ylim)) opz$ylim<-range.ok
      if(is.null(opz$col))  opz$col<-1
      if(is.null(opz$pch))  opz$pch<-1:length(est.psi)
      if(is.null(opz$xlab)) opz$xlab<-"iterations"
      if(is.null(opz$ylab)) opz$ylab<-paste("breakpoint ","(",term,")",sep="")
      if(is.null(opz$type)) opz$type<-"o"
      opz$xaxt<-"n"
      #browser()
      if(is.null(obj$seed)) { #NO boot
        if(all(diff(sapply(obj$psi.history, length)[-1])==0)){ #non-autom (elemento [1] e' NA)
            A<-t(matrix(unlist(obj$psi.history)[-1],nrow=nrow(obj$psi),byrow=FALSE))
            colnames(A)<-rownames(obj$psi)
            opz$x<-0:(nrow(A)-1)
            opz$y<-A[,id.ok]
            par(mfrow=c(1,2))
            do.call(matplot, opz)
            #matplot(0:(nrow(A)-2), A[-1,id.ok],type="o",pch=1:length(est.psi),col=1,
            #      xlab=, ylab=,
            #      ylim=range.ok, xaxt="n",...)
            axis(1,at=0:(nrow(A)-1),cex.axis=.7)
            abline(h=est.psi,lty=3,col=opz$col)
            plot(0:(nrow(A)-1), attr(obj$psi.history,"dev")[-1], ylab="deviance", xlab="iterations", type="o", xaxt="n")
            axis(1,at=0:(nrow(A)-1),cex.axis=.7)
            abline(h = min(attr(obj$psi.history,"dev")),lty=3,col=opz$col)
        } else { #automatic
            psihist<-obj$psi.history[-1]
            id.iter<-rep(1:length(psihist), times=sapply(psihist, length))
            psi.history<-unlist(psihist)
            nomi<-unlist(sapply(psihist, names))
            d<-data.frame(iter=id.iter, psi=psi.history, nomi=nomi)
            #associa i nomi delle componenti di $psi.history (che sono indici 1,2,..) con i nomi della variabile term
            ii<-unique(names(obj$psi.history[[length(obj$psi.history)]])[id.ok])
            if(length(ii)>1) stop("some error in the names?..")            
            with(d[d$nomi==ii,], plot(iter, psi,
                                    xlab=opz$xlab, ylab=opz$ylab,
                                    xaxt="n",...))
            axis(1,at=unique(d$iter),cex.axis=.7)
            #se vuoi proprio associare le stime tra le diverse iterazioni 
            #(per poi unire nel grafico i punti con le linee. Ovviamente alcune linee saranno interrotte)
            #            for(i in 1:length(obj$psi.history)) {
            #                a<-obj$psi.history[[i]]
            
            #                for(j in 1:length(est.psi)){  
            #                    psij<-est.psi[j]
            #a<- ..names match
            #                    r[i,j]<-a[which.min(abs(a-psij))]
            #                    a<-setdiff(a, r[i,j])                    
        }
      } else { #se boot
        par(mfrow=c(1,2))
        plot(obj$psi.history$all.selected.ss, type="b", xlab="bootstrap replicates", 
             ylab="RSS  (selected values)", xaxt="n", pch=20)
        axis(1,at=1:length(obj$psi.history$all.selected.ss),cex.axis=.7)        
        #unicita' delle soluzioni
        if(is.vector(obj$psi.history$all.selected.psi)){
          psi.matr<-m<-matrix(obj$psi.history$all.selected.psi, ncol=1)
        } else {
          psi.matr<-m<-obj$psi.history$all.selected.psi[,id.ok,drop=FALSE]
        }
        
        for(i in 1:nrow(m)) m[i,]<-apply(psi.matr[1:i,,drop=FALSE],2,function(xx)length(unique(xx)))
        m<-t(t(m)+.1*(0:(ncol(m)-1)))
        matplot(1:nrow(m),m, pch=1:ncol(m), type="b", col=1:ncol(m), 
                ylab="no. of distinct solutions",xlab="bootstrap replicates", xaxt="n")
        axis(1,at=1:nrow(m),cex.axis=.7)        
        }
      } #end_fn
