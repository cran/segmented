segreg <- function(formula, data, subset, weights, na.action, family=lm, offset, control=seg.control(), 
    transf=NULL, contrasts=NULL, model=TRUE, x=FALSE, var.psi=TRUE, ...){
  ##### ====================================================================================
  #DA FARE:
  #1) i psi fissi OK prova
  #2) eliminare argomento 'transf'
  #3) se c'e' by consentire un n.psi diverso per categoria di by? ovvero npsi o psi devono essere liste..ok
  #4) la matrice dei contrasti per imporre vincoli alle slope ok 
  #=================================
  #Allora considerando seggrowth() qua ci dovrebbero essere problemi in nel predict..
  
  #=================================
  # `[.withAttributes` <- function(r, i) {
  #   subset <- NextMethod()
  #   attr(subset, "nomeBy") <- attr(r, "nomeBy")
  #   attr(subset, "nomeX") <- attr(r, "nomeX")
  #   attr(subset, "psi") <- attr(r, "psi")
  #   attr(subset, "npsi") <- attr(r, "npsi")
  #   attr(subset, "est") <- attr(r, "est")
  #   attr(subset, "R") <- attr(r, "R")
  #   attr(subset, "fix.psi") <- attr(r, "fix.psi")
  #   attr(subset, "f.x") <- attr(r, "f.x")
  #   attr(subset, "by") <- attr(r, "by")
  #   attr(subset, "levelsBy") <- attr(r, "levelsBy")
  #   subset
  # }
  
    build.all.psi<-function(psi, fixed.psi){
    all.names.psi<-union(names(psi),names(fixed.psi))
    all.psi<-vector("list", length=length(all.names.psi))
    names(all.psi)<- all.names.psi
    for(i in names(all.psi)) {
      if(!is.null(psi[[i]])){
        psi[[i]]<-sort(psi[[i]])
        names(psi[[i]])<-paste("U",1:length(psi[[i]]),".",i,sep="")
      }
      if(!is.null(fixed.psi[[i]])){
        fixed.psi[[i]]<-sort(fixed.psi[[i]])
        names(fixed.psi[[i]])<-	paste("U",1:length(fixed.psi[[i]]),".fixed.",i,sep="")
      }
      all.psi[[i]]<-sort(c(psi[[i]],fixed.psi[[i]]))
    }
    return(all.psi)
  }
  #=================================
    fc<- min(max(abs(control$fc),.8),1)       
    min.step<-control$min.step
    maxit.glm <- control$maxit.glm
    alpha<-control$alpha
    it.max <- old.it.max<- control$it.max
    digits<-control$digits
    toll <- control$toll
    if(toll<0) stop("Negative tolerance ('tol' in seg.control()) is meaningless", call. = FALSE)
    stop.if.error<-control$stop.if.error
    fix.npsi<-fix.npsi<-control$fix.npsi
    if(!is.null(stop.if.error)) {#if the old "stop.if.error" has been used..
      warning(" Argument 'stop.if.error' is working, but will be removed in the next releases. Please use 'fix.npsi' for the future..")
    } else {
      stop.if.error<-fix.npsi
    }
    
    break.boot=control$break.boot
    n.boot<-control$n.boot
    size.boot<-control$size.boot
    gap<-control$gap
    random<-control$random
    pow<-control$pow
    conv.psi<-control$conv.psi
    visual <- control$visual
    visualBoot<-FALSE
    if(visual && n.boot>0) {visual<-FALSE; visualBoot<-TRUE}
    
    
    # if(n.boot>0){
    #   if(!is.null(control$seed)) {
    #     set.seed(control$seed)
    #     employed.Random.seed<-control$seed
    #   } else {
    #     employed.Random.seed<-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
    #     set.seed(employed.Random.seed)
    #   }
    #   if(visual) {visual<-FALSE; visualBoot<-TRUE}# warning("`display' set to FALSE with bootstrap restart", call.=FALSE)}
    #   #        if(!stop.if.error) stop("Bootstrap restart only with a fixed number of breakpoints")
    # }
    last <- control$last
    K<-control$K
    h<-control$h
  
    if (deparse(substitute(family))=="lm" || (is.character(family) && family=="lm")){
      fitter0<-"lm" #get("lm")    
    } else {
      if (is.character(family)) family<-get(family, mode = "function", envir = parent.frame())
      if (is.function(family)) family <- family()
      if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
      }
      fitter0<-"glm"
    }
    
    
    
    s1<-strsplit(as.character(formula)[3],"\\+")[[1]] #separa i termini "additivi"..
    idC<-sapply(sapply(lapply(s1, function(x) grep("seg\\(",x)), function(x) (x>=1)), isTRUE)
    stringa<-s1[idC]  #solo i termini con seg
    
    if(any(sapply(stringa, function(.x) grepl("\\* seg\\(", .x)))) stop("invalid usage of symbol '*' in conjunction with seg()")
    if(any(sapply(stringa, function(.x) grepl("\\:seg\\(", .x)))) stop("invalid usage of symbol ':' in conjunction with seg()")
    
    if(any(sapply(stringa, function(.x) grepl("\\):", .x)))) stop("invalid usage of symbol ':' in conjunction with seg()")
    if(any(sapply(stringa, function(.x) grepl("\\) \\*", .x)))) stop("invalid usage of symbol '*' in conjunction with seg()")

    call <- match.call()
    if (missing(data)) data <- environment(formula)
    tf <- terms(formula, specials = "seg")
    id.ps<-attr(tf,"specials")$seg #posizione nel modelframe; vettore se ci sono piu' termini..include y ma non da interc

    #browser()
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights",  "na.action", "offset"), names(mf), 0L) #"offset",
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    names(mf)[2]<-"formula" #serve se NON hai usato "formula"
    
    #browser()
    
    mf <- eval(mf, parent.frame())
    
    #browser()
    n<-nrow(mf)
    mt <- attr(mf, "terms")
    intercMt<-attr(mt,"intercept")
    interc<-intercMt==1
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) names(Y) <- nm
      }
    if(!is.null(transf)) {
        Y.orig <- Y
        Y <- eval(parse(text=transf), list(y=Y))
        transf.inv<-splinefun(Y, Y.orig, ties=min, method="monoH.FC")
    }
    if(is.null(alpha)) alpha<- max(.05, 1/length(Y))
    if(length(alpha)==1) alpha<-c(alpha, 1-alpha)
    
    #browser()
    .xlivelli<-.getXlevels(mt, mf) 
    weights <- as.vector(model.weights(mf))
    if(!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if(!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
    offs <- as.vector(model.offset(mf))
    
    #browser() #funziona sia nella formula che come argomento?
    
    testo.ps<-names(mf)[id.ps]
    nomiCoefUNPEN<-names(mf)[-c(1,id.ps)]
    # X <- if(!is.empty.model(mt)){
    #     model.matrix(mt, mf, contrasts) 
    #   } else {stop("error in the design matrix")}#matrix(, NROW(Y), 0L)
    # attrContr<-attr(X, "contrasts")
    # n<-nrow(X)
    
    
    #browser()
    #===========================================================================
    #se NON ci sono termini ps
    #===========================================================================
    if(length(testo.ps)<=0) stop("No seg() term. Please, use lm() or glm()")
    #===========================================================================
    #se ci sono termini ps.
    #===========================================================================
    check.estPsi<-function(.x){
      #questo e' utile per verificare che estPsi sia o NA oppure di 0,1 (ma non tutti 0)
      #l'ho fatta perche' se fosse consentito est=c(F,T) (e non FALSE,TRUE) poi T e F venivano presi come nomi di 
      # variabili in predict.segmented
      if(length(.x)==1 && all(is.na(.x))){
        ris<-TRUE
        } else {
          if(!is.numeric(.x)){
            ris<-FALSE
            } else {
              if(length(setdiff(.x,0:1))==0){
                ris<- if(all(.x==0)) FALSE else TRUE
                } else {
                  ris<-FALSE
                }
            }
        }
      ris
    }
    
    
    drop.id<-lambda<-S<-B<-BB<-nomiCoefPEN<-nomiVariabPEN<-NULL
    l<-lapply(testo.ps,function(xx)with(mf,eval(parse(text=xx),envir=data)))
    nomiPS <-unlist(sapply(l,function(xx)attr(xx,"nomeX")))
    psiList <-lapply(l,function(xx)attr(xx,"psi"))
    
    #browser()
    
    #npsiList <-unlist(sapply(l,function(xx) attr(xx,"npsi")))
    npsiList <- lapply(l,function(xx) attr(xx,"npsi"))
    estList <- lapply(l,function(xx) attr(xx,"est"))
    
    #browser()
    
    if(length(setdiff(drop(unlist(sapply(estList, function(.x){if(any(is.na(.x))) 0 else unique(.x)}))), c(0,1)))>0) stop(" 'est' should include 0's and 1's only")
    
    if(length(intersect(nomiCoefUNPEN, nomiPS))>0) stop("any segmented variable included as linear term too?")
    
    #browser()
    
    #if(!all(sapply(estList,check.estPsi))) stop(" 'est' is misspecified in one or more seg() term")

    
    RList <- lapply(l,function(xx) attr(xx,"R"))
    fixpsiList <- lapply(l,function(xx) attr(xx,"fix.psi"))
    fxList <- lapply(l,function(xx) attr(xx,"f.x"))
    byList <- lapply(l,function(xx) attr(xx,"by"))
    nomiBy <- unlist(sapply(l,function(xx) attr(xx,"nomeBy")))
    levelsBy <- lapply(l,function(xx) attr(xx,"levelsBy"))
    
    #browser()
    
    if(all(sapply(levelsBy, is.null)) && (length(npsiList)!=length(nomiPS))) stop(" 'npsi' is not correctly specified")
    
    limZ<-rangeSmooth<-mVariabili<-NULL
		#se ci sono termini ps()+ps(..,by) il nome delle variabili smooth vengono cambiati per aggiungere la variabile by
		nomiPS.orig <- nomiPS
		  
		nomiPS.By <-paste(nomiPS, nomiBy, sep=":")
		nomiPS <- unlist(lapply(nomiPS.By, function(.x) sub(":NULL", "", .x)))

		####se la stessa variabile e' specificata come ps o termine lineare..
		if(length(intersect(nomiPS,nomiCoefUNPEN))>=1) stop("The same variable specified in seg() and as linear term")

		#ATTENZIONE.. se vuoi subito costruire i nomi ps(x), ps(x):z, ecc...usa i:
		nomiPS.ps<- sapply(nomiPS.orig, function(.x)paste("seg(",list(.x),")",sep=""))
		nomiPS.ps<-unlist(lapply(paste(nomiPS.ps, nomiBy, sep=":"), function(.x) sub(":NULL", "", .x)))
		nomiPS.ps.list<-as.list(nomiPS.ps) #serve lista


		for(j in id.ps) mVariabili[length(mVariabili)+1]<-mf[j]
    B<- Bfix <- nomiPS.ps.int.list<-vector(length=length(mVariabili) , "list")
    #BFixed<-BB<-Bderiv
    
    
    #browser()
    
    for(j in 1:length(mVariabili)) {
      if(nomiBy[j]=="NULL"){ # se usuale termine seg()
        nomiBy.j <- NULL
        variabileSmooth<- c(mVariabili[[j]]) #c() converte le matrici in vettori, drop() no..!
        #variabileSmooth<- attr(mVariabili[[j]], "f.x")(variabileSmooth)
        variabileSmooth<- fxList[[j]](variabileSmooth)
        #for(jj in c("nomeX", "psi", "npsi", "f.x", "nomeBy")) attr(variabileSmooth,jj)<-NULL
        B[[j]]<- variabileSmooth
        rangeSmooth[[j]] <- range(variabileSmooth)
        limZ[[j]] <- quantile(variabileSmooth, names=FALSE, probs=c(alpha[1],alpha[2]))
        nomiPS.ps.int.list[[j]]<- nomiPS[j]
        
      } else { #se ci sono termini by
        #browser()
        if(is.null(levelsBy[[j]])){ #se e' vc con variabile continua
          stop(" 'by' in seg(), if provided, should be a factor")
          #B[[j]] <-variabileBy*B[[j]]
          #nomiCoefPEN[[j]]<- sapply(1:ncol(B[[j]]), function(x) gsub(":", paste(".",x, ":", sep="") , nomiPS.ps[j]))
        } else {#se e' VC con variabile categoriale
          nomiBy.j <- nomiBy[j]
          variabileSmooth<- mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE] 
          variabileSmooth<- fxList[[j]](variabileSmooth)
          variabileBy<- mVariabili[[j]][, ncol(mVariabili[[j]]),drop=TRUE]
          M<-model.matrix(~0+factor(variabileBy))
          B[[j]]<- lapply(1:ncol(M), function(.x) M[,.x]*variabileSmooth)
          rangeSmooth[[j]] <- lapply(B[[j]], function(.x) range(.x[.x!=0]))
          limZ[[j]] <- lapply(B[[j]], function(.x) quantile(.x[.x!=0], names=FALSE, probs=c(alpha[1],alpha[2])))
          #browser()
          cond1 <- is.list(psiList[[j]])
          cond2 <- length(names(psiList[[j]]))==length(levelsBy[[j]])
          cond3<- length(setdiff(names(psiList[[j]]),levelsBy[[j]]))==0
          if(cond1&& cond2 && cond3) psiList[[j]]<- psiList[[j]][levelsBy[[j]]]
          nomiPS.ps.list[[j]] <-paste(nomiPS.ps[j], levelsBy[[j]], sep="") #"seg(age):sex1" "seg(age):sex2"
          nomiPS.ps.int.list[[j]]<-gsub("[)]", "", gsub("seg[(]", "", nomiPS.ps.list[[j]])) #age:sexM", "age:sexF"
          #la linea sotto e' se hai diversi breakpoints nei gruppi..
          
        }
      } 
    } #end for(j in 1:length(mVariabili))
    
    #browser()
    #sapply(limZ[[1]], cbind)
    
    repl<-pmax(sapply(B,length)*sapply(B,is.list),1)
    for(i in 1:length(npsiList)){
      if(length(npsiList[[i]])==1) {
        npsiList[[i]] <- rep(npsiList[[i]], repl[i])
        if(!is.list(estList[[i]]) && !is.null(levelsBy[[i]])) estList[[i]] <- rep(estList[i], repl[i])
      }
      if(length(nomiPS.ps.int.list[[i]])!=length(npsiList[[i]])) stop(paste(" 'npsi' (its length) is not correctly specified in the seg term:",i))
      if(!is.null(names(npsiList[[i]]))){
        if(length(setdiff(nomiPS.ps.int.list[[i]], names(npsiList[[i]])))!=0) stop(paste(" 'npsi' (its names) is not correctly specified in the seg term:",i))
      }
    }
    npsiList <- unlist(npsiList)
    
    #browser()
    
    if(!any(sapply(psiList,is.list))) psiList <- rep(psiList, repl)
    if(!any(sapply(estList,is.list))) estList <- rep(estList, repl)
    if(!any(sapply(RList,is.list))) RList <- rep(RList, repl)
    if(!any(sapply(fixpsiList,is.list))) fixpsiList<- rep(fixpsiList, repl)
    
    
    nomiPS.orig <- rep(nomiPS.orig, repl)
    Bfix <- rep(Bfix, repl)
    
    
    
    
    while(any(sapply(B,is.list))){
        id.vc<-which((sapply(B, is.list)))[1]
        nc<-length(B[[id.vc]])
        B<-append(B, B[[id.vc]], after = id.vc-1)
        #for(i in 1:length(B[[id.vc+nc]])) BB<-append(BB, BB[id.vc], id.vc-1)
        B[[id.vc+nc]]<-NULL
        #BB[[id.vc+nc]]<-NULL
        nomiCoefPEN <- append(nomiCoefPEN, nomiCoefPEN[[id.vc]], after = id.vc-1)
        nomiCoefPEN[[id.vc+nc]]<-NULL
        
        psiList <- append(psiList, psiList[[id.vc]], after = id.vc-1)
        psiList[[id.vc+nc]]<-NULL
        
        rangeSmooth <- append(rangeSmooth, rangeSmooth[[id.vc]], after = id.vc-1)
        rangeSmooth[[id.vc+nc]]<-NULL
        
        limZ <- append(limZ, limZ[[id.vc]], after = id.vc-1)
        limZ[[id.vc+nc]]<-NULL
        
        estList <- append(estList, estList[[id.vc]], after = id.vc-1)
        estList[[id.vc+nc]]<-NULL
        
        RList <- append(RList, RList[[id.vc]], after = id.vc-1)
        RList[[id.vc+nc]]<-NULL
        
        fixpsiList <- append(fixpsiList, fixpsiList[[id.vc]], after = id.vc-1)
        fixpsiList[[id.vc+nc]]<-NULL
        
        
        #se la lista contiene solo NULL, non funziona...
        #penMatrixList <- append(penMatrixList, penMatrixList[[id.vc]], after = id.vc-1)
        #penMatrixList[[id.vc+nc]]<-NULL
    }
    
      #if(!all(sapply(estList,check.estPsi))) stop(" 'est' is misspecified in one or more seg() term")
    #browser()
    
    nomiTerminiSEG<-nomiCoefPSI <-NULL
    nomiPS.ps.unlist.seg <- unlist(nomiPS.ps.list)
    nomiPS.ps.unlist <- sub("[)]", "", sub("seg[(]", "",nomiPS.ps.unlist.seg ))
    names(psiList)<- nomiPS.ps.unlist 
    for(i in 1:length(B)) {
        #nomiCoefPSI[[i]]<- paste(paste("psi",1:length(psiList[[i]]), sep=""), nomiPS.ps.unlist[i], sep=".") ##oppure sep=".psi"
        nomiTerminiSEG[[i]]<-rep(nomiPS.ps.unlist[i], length(psiList[[i]]))
    }
    #nomiCoefU<-lapply(nomiCoefPSI, function(.x) sub("psi","U",.x )) 
      #nomiCoefZ<-lapply(nomiCoefPSI, function(.x) sub("psi","Z",.x ))
    nomiSeg<- unique(unlist(nomiTerminiSEG))
      
      
    #browser()
      
      #FINALMENTE (speriamo..:-))
      #nomiCoefZ, nomiCoefpsi, nomiCoefU sono liste con nomi che includono sia le possibili interazioni, sia il n. dei breakpoints
      #Anche nomiTerminiSEG e' della stessa dimensione ma i nomi ignorano il n.dei breakpoints (questa serve per rangeZ)
      
      #nomiSeg
      #npsiList
      #psiList
      #estList
      #RList
      
      # if(is.null(names(estList))) {
      #   names(estList)<-nomiSeg
      # } else {
      #   if(any(sapply(names(estList), function(.x).x==""))) stop(" 'estList' is only partially named. 
      #                                                            Or all or no name allowed.")
      # }
      
    #browser()
      
    if(any(sapply(estList, is.list))) stop(" One or more 'est' components misspecified")
    
    npsiList1<-id.contrR <- rep(NA, length(B))
    psiListE<-psiListQ<-psiList
    for(j in 1:length(B)){
        #K<- npsiList[j]
        K <- if(!is.na(npsiList[nomiSeg[j]])) npsiList[nomiSeg[j]] else npsiList[j]
        npsiList1[j]<- K
        if(any(is.na(psiList[[j]]))){ 
          #if(control$quant) {
            psiListQ[[j]]<- quantile(B[[j]], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)
          #} else {
            psiListE[[j]]<- (min(B[[j]])+ diff(range(B[[j]]))*(1:K)/(K+1))
          #}
        } else {
          K<-npsiList1[j]<-length(psiList[[j]])
        }
        if(!is.null(fixpsiList[[j]])) {
          Bfix[[j]]<- sapply(sort(fixpsiList[[j]]), function(.x) (B[[j]]-.x)*(B[[j]]>.x))
          #colnames(Bfix[[j]])<- paste("U", 1:length(fixpsiList[[j]]),".fixed.",nomiPS.orig[j], sep="")
          colnames(Bfix[[j]])<- paste("U", 1:length(fixpsiList[[j]]),".fixed.",nomiPS.ps.unlist[j], sep="")
          
        }
      #se per qualche termine ci sono le matrici dei vincoli sulle slope
        
        #browser()
        
        j.ok=match(nomiSeg[j], names(RList), nomatch=0)
        j.ok <-if(j.ok>0) j.ok else j 
        if(!any(is.na(RList[[j.ok]]))){
          RList[[j]] <- RList[[j.ok]]
          id.contrR[j] <-TRUE
          } else {
            j.ok=match(nomiSeg[j], names(estList), nomatch=0)
            j.ok <-if(j.ok>0) j.ok else j
            if(!any(is.na(estList[[j.ok]]))){
              if(length(estList[[j.ok]])!=(K+1)) stop(" 'est' is not compatible with 'npsi' ")
              #browser()
              RList[[j]]<-diag(K+1)[,estList[[j.ok]]==1,drop=FALSE]
              id.contrR[j] <-TRUE
            } else {
              RList[[j]]<-diag(K+1)
              id.contrR[j] <-FALSE
            }
          }
    }
    
    #browser()
    
    if(control$quant) {
      psiList<-psiListQ
      initial <- unlist(psiListE)
      PSI1<- matrix(initial, n, length(initial), byrow = TRUE)
      } else {
        psiList<-psiListE
        initial <- unlist(psiListQ)
        PSI1<- matrix(initial, n, length(initial), byrow = TRUE)
      }
    #Quindi PSI1 e' una matrice di valori inziali di psi.. Vanno usati da seg.lm.fit.boot nel caso in cui i primi PSI
    #non portino ad adattare un modello
    
    #NB: Poiche' ora psiList include i veri numeri dei psi, i codici vanno rilanciati
      for(i in 1:length(B)) {
        nomiCoefPSI[[i]]<- paste(paste("psi",1:length(psiList[[i]]), sep=""), nomiPS.ps.unlist[i], sep=".") ##oppure sep=".psi"
        nomiTerminiSEG[[i]]<-rep(nomiPS.ps.unlist[i], length(psiList[[i]]))
      }
      nomiCoefU<-lapply(nomiCoefPSI, function(.x) sub("psi","U",.x )) 
      nomiCoefZ<-lapply(nomiCoefPSI, function(.x) sub("psi","Z",.x ))

      npsii <- sapply(psiList,length)
      id.psi.group <- rep(1:length(psiList), sapply(psiList,length))
      Z<- lapply(1:length(B), function(.x) matrix(B[[.x]], nrow=n, ncol=npsiList1[[.x]]))
      Z<- do.call(cbind,Z)
      colnames(Z) <- unlist(nomiCoefZ)
      
      #nomiPS, nomiPS.By, nomiPS.orig, nomiPS.ps, nomiPS.ps.list, nomiCateg, nomiInterCateg, nomiCoefPEN

      #nomiPS: "x", "z" (vettore)
      #nomiPS.orig: come "nomiPS"
      
      #Se ci sono interazioni (by)
      #nomiPS "x:g"
      #nomiPS.orig: "x", "x", "x".. La stessa variabile ripetuta per il n.dei gruppi
      #
      
      #nomiPS.ps: "seg(x)", "seg(z)" (vettore) [con by: "seg(x):g"]
      #nomiPS.ps.list "seg(x)", "seg(z)" (lista) [lista con by: "seg(x):g1" "seg(x):g2" "seg(x):g3" ..]
      #nomiInterCateg: "x:g1" "x:g2" "x:g3" ..
      
      #========================================================================================================
      X <- if(!is.empty.model(mt)){
        model.matrix(mt, mf, contrasts) 
      } else {stop("error in the design matrix")}#matrix(, NROW(Y), 0L)
      attrContr<-attr(X, "contrasts")
      #n<-nrow(X)
      
      
      X<- X[, !startsWith(colnames(X),"seg("), drop=FALSE]

      idZ <- unlist(tapply(id.psi.group, id.psi.group, function(.x) c(TRUE, rep(FALSE, length(.x)-1))))
      Z.ok<-Z[, idZ, drop=FALSE]
      colnames(Z.ok) <- nomiPS.ps.unlist
      X<-cbind(X, Z.ok) #include anche i termini lineari delle variabili segmented
      #colnames(Z)<- unlist(nomiCoefPEN)
      initial <- unlist(psiList)
      PSI <- matrix(initial, n, length(initial), byrow = TRUE)
      #PSI1 di riserva definito sopra.
      #browser()
      
      #NB la matrie del disegno X include in nomi "seg(x)" e non va bene perche' poi da problemi con i 
      #nomi dei coef dell'oggetto.. Quindi bisogna sostituire questi nomi!!!
      #non serve perche gia' i nomi sono ok..
      #id.segX <-grep( "seg[(]" , colnames(X))
      #colnames(X)[id.segX]<-gsub("[)]", "", gsub("seg[(]", "", colnames(X)[id.segX]))

      if(any(!sapply(Bfix, is.null))){
        X<-cbind(X, do.call(cbind, Bfix))
      }
      
      #browser()
      
      #colnames(X)[unlist(id.psList)] <- nomiPS.orig
      #X[,nomiPS.orig] <- Z[, unique(colnames(Z)), drop=FALSE]
      #nomiCoefPEN include i nomi le interazioni con i livelli (nel caso vc) e anche del numero dei psi
      #[1] "U1.x" "U2.x" "U1.z"
      id.noOW <- if(is.null(weights) && is.null(offs)) TRUE else FALSE 
      if(is.null(weights)) weights<-rep(1,n)
      orig.offs<-offs
      if(is.null(offs)) offs<-rep(0,n)

      #E' realmente necessario assegnare offs e weight con 0 e 1 anche se non servono???
      #browser()
      
      
      limZ <-do.call(cbind, lapply(limZ, function(.x){if(is.list(.x)) do.call(cbind, .x) else cbind(.x)} ))
      #limZ <-matrix(sapply(1:length(npsii), function(.x) rep(limZ[,.x],npsii[.x])), nrow=2, byrow = FALSE)
      limZ <- do.call(cbind, lapply(1:length(npsii), function(.x) matrix(limZ[,.x],nrow=2,ncol=npsii[.x])))
      rangeZ <- do.call(cbind, lapply(1:length(npsii), function(.x) matrix(rangeSmooth[[.x]],nrow=2,ncol=npsii[.x])))
      colnames(rangeZ)<-rep(names(npsii), npsii) #unlist(nomiTerminiSEG)
      
      #browser()
      
      #rangeZ<- matrix(sapply(1:length(npsii), function(.x) rep(rangeSmooth[[.x]],npsii[.x])), nrow=2, byrow = FALSE)
      
      
      invXtX<-Xty<-NULL
      
      
      #browser()
      opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,L0=NULL,visual=visual,it.max=it.max,nomiOK=unlist(nomiCoefU), usesegreg=TRUE,
                fam=family, eta0=NULL, maxit.glm=maxit.glm, id.psi.group=id.psi.group, gap=gap, limZ=limZ, rangeZ=rangeZ,
                conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step, tol.opt=control$tol.opt,
                pow=pow, visualBoot=visualBoot, digits=digits, fc=fc, RList=RList, nomiSeg=nomiSeg, 
                seed=control$seed, min.n=control$min.n, PSI1=PSI1)
      
      #browser()
      
      if(any(id.contrR)){
        if(fitter0=="lm"){
          # for(.i in nomiSeg) { #    #poni min(z)=0, cosi solve() in step.lm.fit non ha problemi.
          #   if(.i %in% colnames(X)) X[,.i]<- X[,.i] - min(X[,.i])
          # }
          
          if(n.boot <= 0) {
          obj <- segConstr.lm.fit(Y, X, Z, PSI, weights, offs, opz)
          } else {
            obj <- segConstr.lm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                                 n.boot = n.boot, size.boot = size.boot, random = random, 
                                 break.boot = break.boot)
  #          seed<- obj$seed
          }
          class0<- "lm"
          if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
        } else {
          if(n.boot<=0){
            obj <-segConstr.glm.fit(Y, X, Z, PSI, weights, offs, opz)
          } else {
            obj <-segConstr.glm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                                   n.boot=n.boot, size.boot=size.boot, random=random, 
                                   break.boot=break.boot)
 #           seed<- obj$seed
          }
          class0<-c("glm","lm")
          }
      } else {
        #browser()

          if(fitter0=="lm"){

          if(n.boot <= 0) {
            obj <- seg.lm.fit(Y, X, Z, PSI, weights, offs, opz)
            } else {
            #browser()
              obj <- seg.lm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                               n.boot = n.boot, size.boot = size.boot, random = random, 
                               break.boot = break.boot)
#              seed<- obj$seed
          }
          class0<-"lm"
          if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
        } else {
          if(n.boot<=0){
            obj <-seg.glm.fit(Y, X, Z, PSI, weights, offs, opz)
          } else {
            obj <-seg.glm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                                n.boot=n.boot, size.boot=size.boot, random=random, 
                                break.boot=break.boot)
   #         seed<- obj$seed
        }
        class0<-c("glm","lm")
        }
      }
      
      #browser()
      
      if(!is.list(obj)){
        warning("Estimation failed. Too many breakpoints? Returning a (g)lm fit..", call. = FALSE)
        
        if(fitter0=="lm"){
          obj0 <- if(id.noOW) lm.fit(x = X, y = Y) else lm.wfit(x = X, y = Y, w = weights, offset = offs)
          class(obj0)<-"lm"
          } else {
            obj0 <- try(suppressWarnings(glm.fit(X, y = Y, offset = offs,
                                               weights = weights, family = opz$fam #control = glm.control(maxit = maxit.glm), 
                                               )), silent = TRUE)
            class(obj0)<-c("glm", "lm")
          }
        return(obj0)
      }
      seed<- obj$seed
      
      if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(invisible(NULL))
      }
      #

      id.psi.group<-obj$id.psi.group
      npsi.groups <- tapply(id.psi.group, id.psi.group, length)
      psi<-obj$psi
      psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
      U<-obj$U
      V<-obj$V
      R=cor(cbind(U,V))
      if(any(abs(R[col(R)>row(R)])>=.99)) stop("Not enought information for the fit: \n too many psi/small sample/replicated data", call.=FALSE)
      rangeZ<-obj$rangeZ
      
      #browser()
      colnames(rangeZ) <- unlist(nomiTerminiSEG)
      it <- obj$it
      epsilon <- obj$epsilon
      id.warn <- obj$id.warn
      k <- length(psi)
      objU <- obj$obj
      #beta.c <- coef(objU)[paste("U", 1:ncol(U), sep = "")]
      beta.c <- coef(objU)[obj$idU]
      
      #browser()
      
      if(any(id.contrR)) {
        beta.c <- lapply(1:length(obj$constr$RList), 
                    function(i) (obj$constr$invA.RList[[i]]%*%beta.c[unlist(obj$constr$nomiUList)==i])[-1])
        beta.c <- unlist(beta.c)
        nomiSeg <- rep(unique(unlist(nomiTerminiSEG)), sapply(obj$constr$nomiUList,function(.x) length(.x)))
        replSeg <- unlist(sapply(obj$constr$nomiUList,function(.x) 1:length(.x)))
        nomiU <- paste(paste("U", replSeg,sep=""), nomiSeg, sep=".")
        nomiVxb <- sub("U","psi", obj$nomiOK)
        X <- obj$X
      } else {
        nomiU <- obj$nomiOK #nomiOK ma puo' essere cambiato se sono eliminati dei psi nella procedura.. 
        nomiVxb <- sub("U","psi", nomiU)
        #In realta' nomiU e nomiVxb gia' ci sono (sarebbero nomiCoefU e nomiCoefPSI), 
        #   pero' se durante la procedura sono stati cambiati perche' alcuni psi sono stati rimossi..bhu??
        #nomiU <-nomiCoefU
        #nomiVxb<- nomiCoefPSI
      }
      
      #browser()
      
      Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
      colnames(U)<- nomiU   #<- nomiOK
      colnames(Vxb)<-nomiVxb #<- sub("U","psi", nomiU)
      se.psi<-rep(NA,k)
      if(fitter0=="lm"){
        objV <- if(id.noOW) lm.fit(x = cbind(X, U, Vxb), y = Y) else lm.wfit(x = cbind(X, U, Vxb), y = Y, w = weights, offset = offs)
        if(any(is.na(objV$coefficients))) {
          warning("some interval with 1 observation. St.errs cannot be computed")
          var.psi=FALSE
        }
        if(var.psi) {
          s2 <- sum(weights*objU$residuals^2)/objV$df.residual
          R <- chol2inv(objV$qr$qr)
          se.psi <- sqrt(diag(R)*s2)[match(nomiVxb, names(coef(objV)),0)]
        }
      } else {
        objV <- try(suppressWarnings(glm.fit(cbind(X, U, Vxb), y = Y, offset = offs,
                                           weights = weights, family = opz$fam, #control = glm.control(maxit = maxit.glm), 
                                           etastart = objU$linear.predictors)), silent = TRUE)
        objV$linear.predictors<-objU$linear.predictors
        objV$deviance<-objU$deviance
        objV$aic<-objU$aic + 2*ncol(PSI) #k
        objV$weights<-objU$weights
        
        #browser()
        
        if (length(offs) && attr(mt, "intercept") > 0L) {
          #se c'e' un offset devi calcolare la null.deviance (come fa glm())
          obj0 <- try(suppressWarnings(glm.fit(X[, "(Intercept)", drop = FALSE], y = Y, offset = offs,
                                               weights = weights, family = opz$fam, #control = glm.control(maxit = maxit.glm), 
                                               etastart = objU$linear.predictors, intercept=TRUE)), silent = TRUE)
          
          # obj0 <- eval(call(if (is.function(method)) "method" else method, 
          #                   x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
          #                   weights = weights, offset = offset, family = family, 
          #                   control = control, intercept = TRUE))
          if (!obj0$converged) warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
          objV$null.deviance <- obj0$deviance
        }
        if(any(is.na(objV$coefficients))) {
          warning("some interval with 1 observation. St.errs cannot be computed")
          var.psi=FALSE
        }
        if(var.psi) {
          R <- chol2inv(objV$qr$qr)
          s2 <- 1
          if(!opz$fam$fam%in%c("poisson","binomial")) s2<- objU$deviance/objV$df.residual
          se.psi <- sqrt(diag(R)*s2)[match(nomiVxb, names(coef(objV)),0)]
        }
      }

      objV$fitted.values <- objU$fitted.values
      objV$residuals <- objU$residuals
      names.coef<-names(objV$coefficients)
      objV$coefficients <- objU$coefficients
      pLin<- ncol(X)
      if(pLin>=1) {
        names(objV$coefficients) <- c(names.coef[1:pLin], c(nomiU, nomiVxb))
      } else {
        names(objV$coefficients) <- c(nomiU, nomiVxb)
      }
      
      #browser()
      ris.psi<-matrix(NA,k,3)
      colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
      rownames(ris.psi) <- nomiVxb
      ris.psi[,2]<-psi
      ris.psi[,3]<- se.psi
      if(stop.if.error)  ris.psi[,1]<-initial
      
      
      all.seg.form<-NULL
      mf1<-mf[1]
      for(i in 2:ncol(mf)) {
        if(i %in% id.ps){
          l<-attributes(mf[[i]])
          if(!is.null(l$by)){
            if(!l$nomeBy%in%names(mf)){
              m<-data.frame(mf[[i]][,1],l$by)
              colnames(m) <- c(l$nomeX, l$nomeBy)
            } else {
              m<-data.frame(mf[[i]][,1])
              colnames(m) <- l$nomeX
            }
            all.seg.form[[length(all.seg.form)+1]]<-as.formula(
              paste("~0+", l$nomeX, "*", l$nomeBy, "-", l$nomeX))
          } else {
            m <-  data.frame(mf[[i]])
            colnames(m) <- l$nomeX
            all.seg.form[[length(all.seg.form)+1]]<- as.formula(paste("~", l$nomeX))
          }
        } else {
          m <-  mf[i]
        }
          mf1<-cbind(mf1, m)
      }
      
      #browser()
      
      #costruisci la formulaLin.. Attenzione non tiene conto di eventuali vincoli sulle pendenze.
      splitFo <- strsplit(as.character(formula),"[+]")
      allX.lin<-paste(c(splitFo[[3]][-grep("seg[(]", splitFo[[3]])], unique(nomiPS.orig)), collapse="+")
      formulaLin <- as.formula(paste(splitFo[[2]], splitFo[[1]], allX.lin))
      
      #browser()
      
      names(all.seg.form)<-nomiPS
      
      #names(mf) <- nomi.mf
      objV$terms <- mt
      if(fitter0=="lm") objV$y<-Y #modificato il 6/5/24.. Non c'e' bisogno.. l'ogg restituit da glm.fit() ha gia' la y
      if(x) objV$x <- X
      objV$contrasts <- attrContr  
      objV$xlevels <- .xlivelli 
      objV$call<-call
      #browser()
      if (model) {
        if(any(!sapply(levelsBy,is.null)) && any(id.contrR)) {
          mf1<- cbind(mf1, Z.ok)
        }
        objV$model <- mf1
      }
      objV$na.action <- attr(mf, "na.action")
      objV$psi <- ris.psi
      objV$id.warn <- id.warn
      objV$it <- it 
      objV$epsilon <- obj$epsilon
      objV$rangeZ <- rangeZ
      objV$constr <- obj$constr
      objV$psi.history <-  psi.values  
      #browser()
      nomiPS.orig<-NULL
      for(i in 1:length(nomiPS)){
        nomiPS.orig[[i]]<-if(is.null(levelsBy[[i]])) nomiPS[i] else paste(nomiPS[i], unlist(levelsBy[[i]]), sep="") 
      }
      nomiPS.orig <- unlist(nomiPS.orig)
        
      #nomiPS.orig<-paste(nomiPS, unlist(levelsBy), sep="")
      
      names(fixpsiList) <- nomiPS.orig
      psi.list<-vector("list", length=length(unique(nomiPS.orig)))
      
      #browser()
      
      names(psi.list)<- nomiPS.orig
      names(psi)<-rep(nomiPS.orig, npsi.groups)
      for(i in names(psi.list)){
        psi.list[[i]]<-psi[names(psi)==i]
      }
      objV$indexU<-build.all.psi(psi.list, fixpsiList)
      #browser()
      objV$nameUV <- list(U = drop(nomiU), V = rownames(ris.psi), Z = nomiPS.orig) #nomiPS.orig??
      objV$nameUV$formulaSeg<- all.seg.form
      objV$nameUV$formulaSegAllTerms<- paste("~", paste(sapply(all.seg.form, function(.x) strsplit(paste(.x), "~"))[2,],collapse="+"))
      objV$formulaLin<- formulaLin
      objV$id.psi.group<- id.psi.group
      objV$psi[,"Initial"]<-NA
      if(n.boot>0) objV$seed <- seed
      #browser()
      
      #il comando structure() l'ho messo perche' avevo bisogno che anche in mancanza di offset, l'oggetto finale restituisse 
      # un oggetto con la componente offset NULL. Cosa che non viene fatta con il semplice comando di sotto.. 
      #objV<- if(id.O) c(objV, offset=offs) else c(objV, offset=NULL)
      #browser()
      objV<- structure(c(objV, list(offset=orig.offs)))
      class(objV)<-c("segmented", class0)
      objV
    
    }
