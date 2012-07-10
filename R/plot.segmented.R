plot.segmented<-function(x, term, add=FALSE, res=FALSE, se=FALSE, show.gap=TRUE,
        link=TRUE, res.col=1, rev.sgn=FALSE, const=0, ... ){
    if (missing(term)) {
        if (length(x$nameUV$Z) > 1) {
            stop("please, specify `term'")
        } else {
            term <- x$nameUV$Z
            }
    } else {
      if(!term%in%x$nameUV$Z) stop("invalid `term'")
      }
    linkinv<-!link
    opz<-list(...)
    cols <- opz$col
    if(length(cols) <= 0) cols <- 1
    lwds <- opz$lwd
    if(length(lwds) <= 0) lwds <- 1
    ltys <- opz$lty
    if (length(ltys) <= 0) ltys <- 1
    cexs <- opz$cex
    if (length(cexs) <= 0) cexs <- 1
    pchs <- opz$pch
    if (length(pchs) <= 0) pchs <- 1
    xlabs<-opz$xlab
    if (length(xlabs) <= 0) xlabs <- term
    ylabs<-opz$ylab
    if (length(ylabs) <= 0) ylabs <- paste("Effect  of ", term, sep=" ")

    a<-intercept(x, term, gap=show.gap)[[1]][,"Est."]
    b<- slope(x,term)[[1]][,"Est."]
    id<-grep(paste("\\.",term,"$",sep=""), rownames(x$psi), value = FALSE)
    est.psi<-x$psi[id,"Est."]
    K<-length(est.psi)
    val<-sort(c(est.psi,x$rangeZ[,term])) #val    avrà K+2 valori
    #a, b avrà K+1 elementi
    a.ok<-c(a[1],a)
    b.ok<-c(b[1],b)
    y.val<-a.ok+b.ok*val + const # a sx di psi
    a.ok1<-c(a,a[length(a)])
    b.ok1<-c(b,b[length(b)])
    #dalla 0.2.9.2 y.val<-y.val1
    y.val<-y.val1<-a.ok1+b.ok1*val + const
    #y.val e y.val1 hanno gli stessi estremi. Tuttavia ci sono difference per i valori in corrispondenza di psi:
    #Per ogni psi y.val sono i valori a sx di psi mentre y.val1 sono quelli a dx..
    s<-1:(length(val)-1)
    xvalues<-x$model[,term]
   	if(rev.sgn) {
		  val<- -val
		  xvalues<- -xvalues
		  }
    m<-cbind(val[s],y.val1[s],val[s+1],y.val[s+1])  #NB ncol=4 e nrow=n.segmenti

    if (inherits(x, what = "glm", which = FALSE) && linkinv){
        fit<- if(res) broken.line(x,term,gap=show.gap,linkinv=linkinv)+resid(x,"response")+const else x$family$linkinv(c(y.val,y.val1))
        xout<-sort(c(seq(val[1],val[length(val)],l=120),val[-c(1,length(val))]))
        l<-approx(as.vector(m[,c(1,3)]), as.vector(m[,c(2,4)]), xout=xout)
        id.group<-cut(l$x, val,FALSE,TRUE)
        yhat<-l$y
        xhat<-l$x
        m[,c(2,4)]<-x$family$linkinv(m[,c(2,4)])
        if(!add) {plot(as.vector(m[,c(1,3)]), as.vector(m[,c(2,4)]),type="n",
            xlab=xlabs, ylab=ylabs, main=opz$main, sub=opz$sub, ylim=range(fit))}
        if(res) points(xvalues, fit,cex=cexs,pch=pchs,col=res.col)
        yhat<-x$family$linkinv(yhat)
        if(length(cols)==1) cols<-rep(cols,max(id.group))
        if(length(lwds)==1) lwds<-rep(lwds,max(id.group))
        if(length(ltys)==1) ltys<-rep(ltys,max(id.group))
        for(i in 1:max(id.group)){
              lines(xhat[id.group==i],yhat[id.group==i],col=cols[i], lwd=lwds[i],lty=ltys[i])
              }
        } else { #end glm
    r<-cbind(val, y.val)
    r1<-cbind(val, y.val1)
    rr<-rbind(r,r1)
    fit<-c(y.val,y.val1)
    if(res) {
        ress<-if (inherits(x, what = "glm", which = FALSE)) residuals(x,"working")*sqrt(x$weights) else resid(x)
        fit<-broken.line(x,term,gap=show.gap,linkinv=linkinv,interc=TRUE)+ress + const
        }
    if(!add) plot(rr, type="n", xlab=xlabs, ylab=ylabs, main=opz$main, sub=opz$sub, ylim=range(fit))
    if(res) points(xvalues, fit,cex=cexs,pch=pchs,col=res.col)
    segments(m[,1], m[,2], m[,3], m[,4], col=cols, lwd=lwds, lty=ltys)
    #for(i in 1:nrow(m)) segments(m[i,1],m[i,2],m[i,3],m[i,4])
    #apply(m, 1, function(z) {segments(z[1], z[2], z[3], z[4])})
    #segments(r[s,1],r[s,2],r[s+1,1],r[s+1,2],lwd=2,col=1:(K+1))
    }
    invisible(NULL)
    }
