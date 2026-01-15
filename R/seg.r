seg <-
function(x, npsi=1, psi=NA, est=NA, R=NA, fixed.psi=NULL, by=NULL, f.x=I){
#------------
    r<-x
    r<- if(!is.null(by)) cbind(r,by) else cbind(r)
    attr(r, "isMatr") <- if(is.matrix(x) && ncol(x)>1) TRUE else FALSE
    nome <- deparse(substitute(x))
    if(is.null(by)) {
      if(is.matrix(x) && ncol(x)>=2) {
        colnames(r) <- if(is.null(colnames(x))) paste(nome, 1:ncol(x), sep="") else colnames(x) 
        } else {
          colnames(r) <-nome
        }
      nome <- colnames(r) 
      #la riga sotto la lascio perche' cosi' come nomeBy restitusce "NULL" piuttostoc he NULL
      attr(r,"nomeBy")<-deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
    } else {
      
      if(is.matrix(by)) {
        if(!is.numeric(by)) stop(" the matrix should be numeric")
        colnames(r)[1]<-nome
        attr(r,"nomeBy")<-paste(colnames(by), collapse=",")
      } else {
        #r<-cbind(x, as.factor(by)) #mettere questo cosi' puo' funzionare anche se by e' "character" o "numeric"?
        colnames(r)<-c(nome, "by") #Perche' non usare attr(r,"nomeBy")
        attr(r,"nomeBy")<-deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
      }
    }
    attr(r,"nomeX")<- nome
    attr(r,"psi")<- psi
    attr(r,"npsi")<- npsi
    attr(r,"est")<- est
    attr(r,"R")<- R
    attr(r,"fix.psi")<- fixed.psi
    attr(r,"f.x")<- f.x
    attr(r, "by")<-by
    attr(r,"levelsBy")<-levels(by)
    #class(r) <- c("withAttributes", class(r))
    r
}


