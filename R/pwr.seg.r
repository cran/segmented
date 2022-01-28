pwr.seg<-function(oseg, pow, n, z="1:n/n", psi, d, s, n.range=c(10,300), X=NULL, break.type=1, alpha=.01, round.n=TRUE, 
                  alternative=c("two.sided","greater","less"), msg=TRUE, ci.pow=0){
  #===Power analysis in segmented regression==
  #Given the input values (z, psi, d..), this function returns n (when pow is provided) or pow (when n is provided)
  #pow, n= the fixed power or sample size. One and only one has to be specified  
  #z= a string indicating the covariate understood to have a segmented effect. Such string should be 
  #   expressed as a function of the quantile function having p as argument.  Namely something like 
  #   qexp(p,..), qbeta(p,..), qunif(p,0,1) or even qunif(p, 1, n) for instance.
  #   It can be a vector meaning the actual covariate, but 'pow' has to be missing. 
  #   Namely if the covariate is supplied (and n is known), it is assumed that only the relevant power can be estimated
  #
  #psi: the breakpoint location (within the covariate range)
  #d= the slope difference
  #s= the response variance
  #n.max= the max sample size to evaluate (ignored if 'n' is provided). However the function can estimate sample sizes 
  #   larger than n.max
  #X: the design matrix including additional linear variables
  #alpha= the significance level
  #round.n= if TRUE the 'estimated' sample size is rounded 
  #authors: Nicoletta D'Angelo and Vito Muggeo
  #---------------
  #Examples:
  ##  pwr.seg(n=100,psi=.5, d=2, s=1)
  ##  pwr.seg(pow=.23,psi=5, d=.5, s=5)
  ##  pwr.seg(pow=.23,z="1:n/n",psi=.5, d=1.5, s=1)
  pwr<-function(x, psi, d, s, X=NULL, alpha=.01, alt=c("two.sided","greater","less"), change=1){
    n<-length(x)
    v<-quantile(x, c(.02,.98), names=FALSE)
    if(psi<=v[1] || psi>=v[2]) stop("psi outside the covariate range", call.=FALSE)
    values <- seq(min(x), max(x), length = 20)
    n1 <- length(values)
    PSI <- matrix(values, nrow = n, ncol = n1, byrow = TRUE)
    X1 <- matrix(x, nrow = n, ncol = n1, byrow = FALSE)
    
    if(change==1) {
      X2<- (X1-PSI)*(X1>PSI)
      u.psi<-((x-psi)*(x>psi))
    } else {
      X2<-1*(X1>PSI)
      u.psi<- 1*(x>psi)
    }
    pmaxMedio <- rowMeans(X2)
    #pmaxMedio <- ((x-psi)*(x>psi))
    if(is.null(X)) X<-cbind(1,x)
    ImH <- -X %*% solve(crossprod(X),t(X))
    diag(ImH)<- 1+ diag(ImH)
    th <- d * pmaxMedio %*% ImH %*% u.psi
    #se <- s * drop(sqrt(pmaxMedio %*% ImH %*% pmaxMedio))
    se <- s* sqrt(rowSums(pmaxMedio %*% ImH * pmaxMedio))
    #
    #browser()
    pw<-switch(alternative, 
               "greater" =  pnorm(qnorm(alpha, sd = se, lower.tail = FALSE),
                            mean = th, sd =  se, lower.tail = FALSE),
               "less"    =  pnorm(qnorm(alpha, sd = se, lower.tail = TRUE),
                            mean = th, sd =  se, lower.tail = TRUE),
               "two.sided"= pnorm(qnorm(alpha/2, sd = se, lower.tail = TRUE),
                               mean = th, sd =  se, lower.tail = TRUE) + 
                          pnorm(qnorm(alpha/2, sd = se, lower.tail = FALSE),
                            mean = th, sd =  se, lower.tail = FALSE)
        )
    pw
  }
  #========================================
  #browser()
  x<-z
  alternative <- match.arg(alternative)
  if(!(break.type %in% 1:2)) stop(" 'break.type' should be 1 or 2")
  if(missing(oseg)){
      if(!is.character(x) && !missing(pow)) stop(" if the covariate is provided, 'pow' has to be missing ")
      if(!is.character(x) && !missing(n))   stop(" if the covariate is provided, 'n' has to be missing ") 
    
      if(missing(pow)) {
        if(is.character(x)){
          if(missing(n)) stop(" 'n' or 'pow' have to be provided")
          p<-seq(.001,.999,l=n)
          x<-eval(parse(text= x))
        } else {
          n<-length(x)
        }
        pow<-pwr(x=x, psi=psi,d=d, s=s,X=X, alpha=alpha, alt=alternative, change=break.type)
        if(msg) cat("Est. power:", paste(round(pow,3)), "\n")
          else return(pow)
      } else {
        if(!missing(n)) stop(" Only 'n' *or* 'pow' can be provided")
        if(pow<=.0001 || pow>=.9999) stop(" 'pow' should be in (.0001, .9999)")
        n.values<-round(seq(min(n.range), max(n.range), l=50))
        K<-length(n.values)
        pp<-rep(NA, K)
        x0<-x
        for(i in 1:K) {
          n<-n.values[i]
          p<-seq(.001,.999,l=n)
          x<-eval(parse(text= x0))
          pp[i]<-pwr(x=x, psi=psi, d=d, s=s, X=X, alpha=alpha, alt=alternative, change=break.type)
        }
        #browser()
        if(length(unique(pp))<=3) {
          cat(stop("Too few distinct values (", unique(pp), ") in the computed power(s)\n", call.=FALSE))
        }
        r<-cbind(n.values, pp)
        f<-splinefun(r[,2], r[,1], method="monoH.FC") #method="hyman"
        a<-f(pow)
        if(round.n) a<-round(a)
        if(msg) cat("Est. sample size:", paste(a), "\n")
          else return(a)
      }
  } else { #se c'e' l'ogg "segmented"
      if(length(oseg$nameUV$V)>1) stop("only models with just 1 breakpoint allowed")
      #introdurre argomento id.psi??? oseg$nameUV$V[id.psi]
      #browser()  
      x<-oseg$model[,oseg$nameUV$Z]
      n<-length(x)
      psi<-oseg$psi[oseg$nameUV$V,2]
      d<-coef(oseg)[oseg$nameUV$U]
      s<-summary(oseg)$sigma
      X<-model.matrix(oseg)
      X<-X[, setdiff(colnames(X),c(oseg$nameUV$U, oseg$nameUV$V)),drop=FALSE]
      pow<-pwr(x=x, psi=psi, d=d, s=s, X=X, alpha=alpha, alt=alternative, change=break.type)
      
      if(ci.pow>0){
        V<-vcov(oseg)[c(oseg$nameUV$U,oseg$nameUV$V),c(oseg$nameUV$U,oseg$nameUV$V)]
        powi<-rep(NA,ci.pow)
        for(i in 1:ci.pow){
          r<-MASS::mvrnorm(1, c(d,psi), V)
          di<-r[1]
          psii<-r[2]
          powi[i]<-pwr(x=x, psi=psii, d=di, s=s, X=X, alpha=alpha, alt=alternative, change=break.type)
        }
        ci<-quantile(powi,c(.025,.975),names=FALSE)
        m<- paste("Est. power for the current fit:", paste(round(pow[1],3)), " (", paste(round(ci,3), collapse = ", "),")")
        pow<-c(pow,ci)
      } else {
        m <- paste("Est. power for the current fit:", paste(round(pow[1],3)))
      }
      if(msg) cat(m, "\n")
      #if(msg) cat("Est. power for the current fit:", paste(round(pow,3)), 
       #           " (", paste(round(ci,3), collapse = ", "),")", "\n")
       else return(pow)
  }
}

