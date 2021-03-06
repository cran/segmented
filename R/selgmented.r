selgmented<-function(olm, seg.Z, alpha=0.05, type=c("score" , "davies", "bic"), control=seg.control(),
                          return.fit=TRUE, bonferroni=FALSE, Kmax=2, msg=TRUE){
  #Selecting number of breakpoints in segmented regression (via the R package segmented)
  #Author: vito.muggeo@unipa.it
  if(missing(seg.Z)){
    nomeX <- all.vars(formula(olm))[2]
    if(length(nomeX)>1 || any(is.na(nomeX))) stop("I cannot determine the segmented variable")
    seg.Z<- as.formula(paste("~", nomeX ))
  } else {
    if(length(all.vars(seg.Z))>1) stop("Multiple variables are not allowed in seg.Z")
  }
  type<-match.arg(type)
  if(type!="bic" && Kmax!=2) stop("Kmax>2 is not (yet) allowed with hypothesis testing procedures", call.=FALSE)
  if(type=="bic"){
    npsi<-1:Kmax
    ris<-vector("list", length(npsi))  
    bic.values<-rep(NA, length(npsi))
    for(i in npsi){
      ris[[i]]<-suppressWarnings(try(segmented(olm, seg.Z, npsi=i, control=control, silent=TRUE)))
      if(!inherits(ris[[i]], "segmented")) ris[[i]]<-suppressWarnings(try(segmented(olm, seg.Z, npsi=i, control=control, silent=TRUE))) 
      if(inherits(ris[[i]], "segmented")) bic.values[i]<- BIC(ris[[i]]) #-2*logLik(ris[[i]]))+ edf*log(n)*Cn
      #if(inherits(ris[[i]], "segmented")) {
      #  osum<-summary(ris[[i]])
      #  edf<-osum$df[1]
      #  n<-sum(osum$df[1:2])
      #  sigma2<-sum(ris[[i]]$residuals^2)/n
      #  bic.values[i]<- log(sigma2)+ edf * (log(n)/n)*1 
      #}
    }
    bic.values<-c(BIC(olm), bic.values)
    names(bic.values)<-c("0", npsi)
    n.psi.ok<- c(0,npsi)[which.min(bic.values)]
    r<-list(bic.values=bic.values, n.psi=n.psi.ok)
    if(!return.fit) {
      return(r) #return(c(n.psi.ok, NA))
    }
    
    if(msg){
      cat("BIC to detect no. of breakpoints\n")
      cat("BIC values:\n")
      print(bic.values)
      cat(paste("No. of selected breakpoints: ", n.psi.ok, " \n"))
    }
    ris[[n.psi.ok]]$selection.psi <- bic.values
    return(ris[[n.psi.ok]])
  }
  alpha.adj<-alpha/Kmax
  p1<- if(type=="score")  pscore.test(olm, seg.Z, n.break=2)$p.value else davies.test(olm)$p.value
  p1.label<-"p-value '0 vs 2' "
  if(p1>alpha.adj){
    p2.label<-"p-value '0 vs 1' "
    p2<- if(type=="score") pscore.test(olm, seg.Z, n.break=1)$p.value else p1 #davies.test(olm)$p.value
    if(!bonferroni) alpha.adj<- alpha
    if(p2>alpha.adj) {
      out<-olm
    } else {
      out<-segmented(olm, seg.Z, npsi=1, control=control)
    }
  } else {
    p2.label<-"p-value '1 vs 2' "
    #################
    olm<-update(olm, data=model.frame(olm)) #questo e' necessario per far funzionare davies.test() sotto..
    ################
    o1<-segmented(olm, seg.Z, npsi=1, control=control)
    if(type=="score") {
      p2<-pscore.test(o1, seg.Z, more.break=TRUE)$p.value
    } else {
      #KK<-new.env()
      #olm1<-update(olm, data=model.frame(o1))
      #o1<-  update(o1, obj=olm1)
      p2<-  davies.test(o1, seg.Z)$p.value
    }
    if(!bonferroni) alpha.adj<-alpha 
    if(p2>alpha.adj) {
      o1<-segmented(olm, seg.Z, npsi=1, control=control)
      #cat("One breakpoint detected\n")
      out<-o1
    } else {
      o2<-segmented(olm, seg.Z, npsi=2, control=control)
      #cat("Two breakpoint detected\n")
      out<-o2
    }
  }
  n.psi.ok<-length(out$psi[,2])
  x2<- -2*sum(log(c(p1,p2)))
  p<-1-pchisq(x2, df=2*2)
  r<-list(pvalues=c(p1=p1, p2=p2, p=p), npsi=n.psi.ok)
  attr(r, "label")<- p2.label
  if(!return.fit) {
    return(r)
  }
  if(msg){
    cat("Hypothesis testing to detect no. of breakpoints\n")
    cat(paste("statistic:", type,"  level:", alpha, "  Bonferroni correction:", bonferroni, "\n"))
    cat(paste(p1.label, "= ", format.pval(p1,4), "   ", p2.label, "= ", format.pval(p2,4) ,
            " \nOverall p-value = ", format.pval(p,4),"\n",sep=""))
    cat(paste("No. of selected breakpoints: ", n.psi.ok, "\n"))
  }
  out$selection.psi<-r
  return(out)
}
