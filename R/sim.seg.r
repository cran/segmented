
sim.seg<-function(N, ni, fixEff=c(2,-.2,.3,.4),
            S=diag(c(.05,.05,.05,.07)^2), seg.Z=NULL,
            psii, fam=c("gauss","pois","bin"), 
            z.psi, x.diff,  #x.diff=list(.1,rbinom(N,size=1,prob=.5) )
            sd.eps=.01, seed, display=FALSE, eta0=FALSE){
#VERSIONE: 08/10/2009 - 4/10/2010 - 30/4/18 - 25/06/19 
#Funzione per simulare dati brokenline con effetti casuali
#--arguments
#N: number of subjects
#n: the cluster size. If scalar, 'n' measurements per subject are generated, otherwise the N actual cluster sizes are computed as
#     sample(min(n):max(n), N, replace=TRUE).  
#fixEff: the fixed effects (a,b,c,psi) into a+b*x+c*(x-psi,0) (che rappresentano anche le medie dei random effects). 
#   3 or 4 values can be provided, the fourth representing the breakpoint. If psii is provided, fixEff should have 3 values.
#S: la matrice di var per i tre (o quattro) effetti casuali a,b,c (ed eventualmente psi)
#x: the segmented covariate for each subject. It is a string (possibly as a function of n) such as the actual values  "c(2,5,6.7)" or 
#    "runif(n,0,1)" or "1:n". Note the actual x-values are computed as "rep(x, length.out=n)", therefore shorter or longer vectors are (tacitly)  #     accepted. if NULL, it is taken sort(runif(n,0,1))
#psii: if provided, it represents the vector of the actual breakpoint values. If missing, the breakpoints are generated via a normal distribution

  
  
  #psi: a vector with three components: the fixed effect for the psi parameter, a string to simulate the ranef and the link function
#     If length(fixEff)=4, this argument is ignored.. Of course the physical meaning depends on the link function and on z.psi
#     For instance, if z.psi is provided, psi[1] is the intercept, namely the breakpoint when z=0
#z.psi, x.diff: Eventuali liste di due componenti: il coefficiente e la variabile. Ad es.,
#   x.diff=list(.1,rbinom(n,size=1,prob=.5)) significa che la DiffSlope dipende da .1*rbinom(n,size=1,prob=.5)
#   Se mancante nessuna esplicativa viene considerata
#sd.e: la dev.st della risposta
#display: se TRUE i dati vengono anche disegnati
#Returns a list with four components:
#   i)data=the simulated data, including the variables z.psi and x.diff with names "X.diff" and "Z.psi"
#   ii)coef=matrix of coeffs (random + fixed) "b0i","b1i","di","kappa0i"
#   iii)randeff= matrix of randeff "b0i" "b1i" "di"  "ki"
#   iv)fixeff: the fixed effects
#   v)psii: the actual breakpoints, possibly plogis() of (kappa0+k_i+kappa1*Z_i)
#   vi)var: the random effect covariance matrix
#author: vito muggeo
  get_rand_state <- function() {
    # Using `get0()` here to have `NULL` output in case object doesn't exist.
    # Also using `inherits = FALSE` to get value exactly from global environment
    # and not from one of its parent.
    get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  
  set_rand_state <- function(state) {
    # Assigning `NULL` state might lead to unwanted consequences
    if (!is.null(state)) {
      assign(".Random.seed", state,  envir = .GlobalEnv,  inherits = FALSE)
    }
  }
  
  
  #old_state <- get_rand_state()
  #on.exit(set_rand_state(old_state))
  
  #------------
  y<-id<-NULL
  #require(MASS) #per mvrnorm()
  fam<-match.arg(fam)
  f.sim<- function(etai,fam){
    ni=length(etai)
    switch(fam,
           gauss=rnorm(ni, etai, sd.eps),
           pois=rpois(ni,exp(etai)),
           bin=rbinom(ni,1,plogis(etai))
    )
  }
  #browser()
  #psi<-psii
  n<-ni
  if(length(ni)<=1 && ni==1){
    idFix=TRUE
    S=diag(rep(0,4))
  } else {
    idFix=FALSE
  }
  logit<-function(x)log(x/(1-x))
  if(length(fixEff)!=ncol(S)) stop("length(fixEff) does not match dim(S)")
  if(!missing(seed)) set.seed(seed)
  RANEF<-mvrnorm(n = N, mu=rep(0,ncol(S)), Sigma=S, empirical=FALSE)
  COEF<-t(t(RANEF) + fixEff)
  
  #browser()
  
  if(!missing(psii)){
    if(!(is.vector(psii) && length(psii)==N)) stop("if psii is provided, it should be length(psii)=N")
    if(length(fixEff)>3) stop("if 'psii' is provided, 'fixEff' should have 3 values")
    link<-psi.link<-"identity"
    if(!link%in%c("identity","logit")) stop("invalid psi.link")
    COEF<-cbind(COEF,1)
    RANEF<-cbind(RANEF,1)
    #kappa0<-eval(parse(text=psi[1]))
    #ki<-RANEF[,4]<-eval(parse(text=psi[2]))
    kappa0<-mean(psii)
    ki<- psii-kappa0
    kappa0i<-COEF[,4]<-kappa0+ki
    } else {
      link<-"identity"
      kappa0<-fixEff[4]
      kappa0i<-COEF[,4]
    }
  colnames(COEF)<-c("beta0i","beta1i","deltai","kappa0i")
  colnames(RANEF)<-c("b0i","b1i","di","ki")
  #esplicative nella diffSlope e psi
  kappa1<-delta1<-0
  #xx.diff<-zz.psi<-rep(0,N)
  id.x.diff<-id.z.psi<-FALSE
  #browser()
  eta.diff<-0
  if(!missing(x.diff)){
    if(is.numeric(x.diff)&& length(x.diff)==1) x.diff<-list(x.diff[1],rep(c(0, 1), length = N) )
    if(!(is.list(x.diff) && length(x.diff)==2))stop(" 'x.diff' has to be a list ")
    id.x.diff<-TRUE
    id.coef <- which.min(sapply(x.diff, length))
    delta1<-x.diff[[id.coef]]
    xx.diff<- x.diff[[-id.coef]]
    if(length(delta1)>1){
      if(!is.matrix(xx.diff)) stop("xx.diff should be a matrix")
      eta.diff<-drop(xx.diff%*%delta1)
    } else {
      if(!(is.vector(xx.diff) && length(xx.diff)==N)) stop("x.diff must be a N-vector")
      eta.diff<-delta1*xx.diff
    }
    
  }
  eta.psi<-0
  #browser()
  if(!missing(z.psi)){
    if(is.numeric(z.psi)&& length(z.psi)==1) z.psi<-list(z.psi[1], rep(c(0,1), c(ceiling(N/2), floor(N/2))))
    if(!(is.list(z.psi) && length(z.psi)==2)) stop(" 'z.psi' has to be a list ")
    id.z.psi<-TRUE
    id.coef <- which.min(sapply(z.psi, length))
    kappa1<-z.psi[[id.coef]]
    zz.psi<-z.psi[[-id.coef]] #eval(parse(text=z.psi[[2]]))
    if(length(kappa1)>1){
      if(!is.matrix(zz.psi)) stop("should be a matrix")
      eta.psi<-drop(zz.psi%*%kappa1)
    } else {
      #if(length(zz.psi)!=N) stop("z.psi must be a N-vector")
      if(!(is.vector(zz.psi) && length(zz.psi)==N)) stop("zz.psi must be a N-vector")
      eta.psi<-kappa1*zz.psi
    }
  }
  if(link=="identity"){
    psi.ok<-eta.psi + kappa0i
    psiFix<- kappa0
    } else {
      psi.ok<-plogis(eta.psi + kappa0i)
      psiFix<- plogis(kappa0)
    }

  if(idFix){
    if(is.null(seg.Z)) {
      xi<- if(is.character(seg.Z)) eval(parse(text=seg.Z)) else sort(runif(N,0,1))
    } else {
      if((is.vector(seg.Z)&&length(seg.Z)==N)) xi<-seg.Z 
    }
    Xi<-cbind(1,xi,pmax(xi-psi.ok,0))
    etai<-drop(Xi%*%COEF[1,1:3])+ eta.diff*pmax(xi-psi.ok,0)
    yi<- f.sim(etai,fam) #etai + rnorm(ni[i],0,sd.e)
    Ris<-data.frame(y=yi,x=xi, eta0=etai)
    if(id.x.diff) Ris$x.diff<- xx.diff
    if(id.z.psi) Ris$z.psi<- zz.psi 
    
  } else {
    if(length(unique(n))==1) ni <- rep(unique(n),N) else ni<-sample(min(n):max(n), size=N, replace=TRUE)
    eta.diff<-rep(eta.diff, length=N)
    Eta<-X<-Y<-list()
    for(i in 1:N){
      n<-ni[i]
      xi<- if(is.null(seg.Z)) sort(runif(n,0,1)) else  eval(parse(text=seg.Z))
      xi<-rep(xi, length=n)
      xi<-sort(xi)
      Xi<-cbind(1,xi,pmax(xi-psi.ok[i],0))
      etai<-drop(Xi%*%COEF[i,1:3])+ eta.diff[i]*pmax(xi-psi.ok[i],0) #delta1*xx.diff[i]*pmax(xi-psi.ok[i],0)
      yi<- f.sim(etai,fam) #etai + rnorm(ni[i],0,sd.e)
      Eta[[length(Eta)+1]]<-etai
      X[[length(X)+1]]<-xi
      Y[[length(Y)+1]]<-yi
    }
    Ris<-data.frame(cbind(id=rep(1:N,times=ni), y=unlist(Y),x=unlist(X)), eta0=unlist(Eta))
    if(id.x.diff) {
      if(is.matrix(xx.diff)) {
        xdiffD<-do.call(rbind, lapply(1:nrow(xx.diff), function(.x) matrix(xx.diff[.x,], nrow=ni[.x], byrow=TRUE,ncol=ncol(xx.diff))))
        colnames(xdiffD)<-paste("Xdif", 1:ncol(xdiffD),sep="")
        Ris=cbind(Ris, xdiffD)
      } else {
        Ris$Xdif<-rep(xx.diff,times=ni)
      }
    }
    if(id.z.psi) {
      if(is.matrix(zz.psi)) {
        zpsiD<-do.call(rbind, lapply(1:nrow(zz.psi), function(.x) matrix(zz.psi[.x,], nrow=ni[.x], byrow=TRUE,ncol=ncol(zz.psi))))
        colnames(zpsiD)<-paste("Zpsi", 1:ncol(zpsiD),sep="")
        Ris<-cbind(Ris, zpsiD)
      } else {
        Ris$Zpsi<- rep(zz.psi,times=ni)
      }
    }
  }

  #browser()
  
  if(display) {
    xx<-seq(min(Ris$x),max(Ris$x),length=100)
    #browser()
    f1<- function(.x,fam){
      switch(fam,
             gauss=.x,
             pois=exp(.x),
             bin=plogis(.x)
      )
    }
    
    segg<-function(.x,b) b[1]+b[2]*.x+b[3]*pmax(.x-b[4],0)
    if(idFix){
      plot(y~x, data=Ris, pch=20, ylim=range(y,fixEff[1]+fixEff[2]*kappa0))
      if(!id.z.psi && !id.x.diff){
        lines(f1(segg(sort(c(range(Ris$x), psiFix)), c(fixEff[1:3],psiFix)),fam) ~sort(c(range(Ris$x), psiFix)) , data=Ris, lwd=4)
        rug(psiFix, lwd=4)
      }
      
    } else {
      plot(y~x, data=Ris, type="n")
      for(i in 1:N) {
        Risi<-subset(Ris, id==i)
        with(Risi, points(x,y,col=i,pch=20))
        xxi<-sort(c(range(Risi$x), psi.ok[i]))
        XXi<-cbind(1,xxi,pmax(xxi-psi.ok[i],0))
        etai<-drop(XXi%*%COEF[i,1:3])+ eta.diff[i]*pmax(xxi-psi.ok[i],0)
        lines(xxi, f1(etai,fam), col=i) 
        rug(psi.ok[i], col=i)
        #questo ? per aggiungere la relazione segmented di effetti fissi..
        #psi.i<-psi.ok[i]
                  #segFix.i<- fixEff[1]+fixEff[2]*xx+(fixEff[3]+delta1*xx.diff[i])*pmax(xx-psi.i,0)
                  #points(xx, segFix.i, type="l",lwd=2)
      }
      if(!id.z.psi && !id.x.diff){
        lines(f1(segg(sort(c(range(Risi$x), psiFix)), c(fixEff[1:3],psiFix)),fam) ~sort(c(range(Risi$x), psiFix)) , data=Ris, lwd=4)
        rug(psiFix, lwd=4)
      }
    }
  }
  #browser()
  
  FixEff<-  c(fixEff[1:3], if(id.x.diff) delta1 else NULL, kappa0, if(id.z.psi) kappa1 else NULL)
  nnomi<- c("beta0","beta1","delta0", if(id.x.diff) paste("delta", 1:length(delta1),sep="") else NULL, "kappa0", if(id.z.psi) paste("kappa", 1:length(kappa1),sep="") else NULL)
  names(FixEff)<-nnomi #[1:3]<-c("beta0","beta1","delta0")
  #call<-match.call()
  sd.eps<-ifelse(fam=="gauss", sd.eps, 1)
  if(!eta0) Ris$eta0<-NULL
  Ris<-list(data=Ris,fixef=FixEff,psi.i=psi.ok,sd=sd.eps) #,call=call) #,perc.psi=perc.psi)
  if(!idFix){
    perc.psi<-tapply((Ris$data$x<rep(psi.ok, ni)), Ris$data$id, mean)
    Ris$ni<-ni
    Ris$perc.psi<-perc.psi
    Ris$ranef=RANEF
    Ris$sd=list(sqrt(diag(S)),sd.eps)
    #coef=COEF,
  }
  #browser()
  
  return(Ris)
}


