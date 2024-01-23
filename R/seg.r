seg <-
function(x, npsi=1, psi=NA, est=NA, R=NA, fixed.psi=NULL, by=NULL, f.x=I){
#------------
    #browser()
    r<-x
    nome <- deparse(substitute(x))
    r<- if(!is.null(by)) cbind(r,by) else cbind(r)  
    colnames(r)<- if(!is.null(by)) c(nome, "by") else nome
    attr(r,"nomeX")<- nome
    attr(r,"psi")<- psi
    attr(r,"npsi")<- npsi
    attr(r,"est")<- est
    attr(r,"R")<- R
    attr(r,"fix.psi")<- fixed.psi
    attr(r,"f.x")<- f.x
    attr(r, "by")<-by
    attr(r,"levelsBy")<-levels(by)
    attr(r,"nomeBy")<-deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
    r
}


#segmented.formula(y~seg(z, npsi=2, by=sex)+sex+seg(age))
#segmented.formula(y~seg(age)+ seg(z, npsi=2, by=sex)+sex)
#segmented.formula(y~seg(age, by=sex)+ seg(z, npsi=2, by=sex)+sex)
