#funziona fino a un certo punto.. 
#E comunque devi aggiunegere il caso di interazioni...

#byList nomiBy sono definiti.. 

#Vedi C:\dati\lavori\jumpoint\fasola\funzioniSalvo\nuove\perpacchetto

#stepreg(y~seg(tt, by=x1, npsi=2)+seg(tt, by=cbind(x2,x3))+seg(tt))
stepreg <- function(formula, data, subset, weights, na.action, family=lm, control=seg.control(), 
    transf=NULL, contrasts=NULL, model=TRUE, x=FALSE, var.psi=FALSE, ...){
  ##### ====================================================================================
  #DA FARE:
  #1) i psi fissi OK prova
  #2) eliminare argomento 'transf'
  #3) se c'e' by consentire un n.psi diverso per categoria di by? ovvero npsi o psi devono essere liste..ok
  #4) la matrice dei contrasti per imporre vincoli alle slope ok 
  #=================================
  #Allora considerando seggrowth() qua ci dovrebbero essere problemi in nel predict..
  
  # ---------
  #a differenza delle funzioni seg*, nei modelli step*, il modello finale da cui estrarre i fitted viene adattato 
  #NON nelle funzioni step.lm.fit, ma in stepmented.(g)lm.. Quindi facciamo lo stesso in stepreg() (anche se non dovrebbe essere cosi'--)
  step.lm.fitSC.boot<-step.lm.fitSC<-step.glm.fitSC.boot<-step.glm.fitSC<-NULL
  stepConstr.lm.fit.boot<- stepConstr.lm.fit<- stepConstr.glm.fit.boot<-stepConstr.glm.fit<-NULL
  mylm<-function(x,y,w=1,offs=0){
    x1<-x*sqrt(w)
    y<-y-offs
    y1<-y*sqrt(w)
    XtX <- crossprod(x1)
    b<-drop(solve(XtX,crossprod(x1,y1)))
    fit<-drop(tcrossprod(x,t(b)))
    r<-y-fit
    o<-list(coefficients=b,fitted.values=fit,residuals=r, df.residual=length(y)-length(b), XtX=XtX, w=w)
    o
  }
  #-----------
  #=================================
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
    #===============================
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
    #===============================
    
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
    
    break.boot=control$break.boot + 2
    n.boot<-control$n.boot
    size.boot<-control$size.boot
    gap<-control$gap
    random<-control$random
    pow<-control$pow
    conv.psi<-control$conv.psi
    display <- control$visual
    #visualBoot<-FALSE
    #if(visual && n.boot>0) {visual<-FALSE; visualBoot<-TRUE}
    agg<- 1-control$fc
    
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

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights",  "na.action"), names(mf), 0L) #"offset",
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    names(mf)[2]<-"formula" #serve se NON hai usato "formula"
    
    mf <- eval(mf, parent.frame())
    
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
    
    testo.ps <-names(mf)[id.ps]
    nomiCoefUNPEN<-names(mf)[-c(1,id.ps)]
    X <- if(!is.empty.model(mt)){
        model.matrix(mt, mf, contrasts) 
      } else {stop("error in the design matrix")}#matrix(, NROW(Y), 0L)
    attrContr<-attr(X, "contrasts")
    n<-nrow(X)
    
    #browser()
    #===========================================================================
    #se NON ci sono termini ps
    if(length(id.ps)<=0) stop("No seg() term. Please, use lm() or glm()")
    #se ci sono termini ps.
    
    drop.id<-lambda<-S<-B<-BB<-nomiCoefPEN<-nomiVariabPEN<-NULL
    l<-lapply(testo.ps,function(xx)with(mf,eval(parse(text=xx),envir=data)))
    nomiPS <- drop(unlist(sapply(l,function(xx)attr(xx,"nomeX"))))
    psiList <-lapply(l,function(xx)attr(xx,"psi"))
    npsiList <- lapply(l,function(xx) attr(xx,"npsi"))
    estList <- lapply(l,function(xx) attr(xx,"est"))
    if(length(setdiff(drop(unlist(sapply(estList, function(.x){if(any(is.na(.x))) 0 else unique(.x)}))), c(0,1)))>0) stop(" 'est' should include 0's and 1's only")
    
    #if(length(intersect(nomiCoefUNPEN, nomiPS))>0) stop("any segmented variable included as linear term too?")
    
    RList <- lapply(l,function(xx) attr(xx,"R"))
    fixpsiList <- lapply(l,function(xx) attr(xx,"fix.psi"))
    fxList <- lapply(l,function(xx) attr(xx,"f.x"))
    byList <- lapply(l,function(xx) attr(xx,"by"))
    nomiBy <- unlist(sapply(l,function(xx) attr(xx,"nomeBy")))
    levelsBy <- lapply(l,function(xx) attr(xx,"levelsBy"))
    names(byList) <- nomiBy
    
    #browser()
    ## questo e' un tentativo per tenere conto di seg(X, npsi=c(2,1,2...))
    if(length(nomiPS)>1 && length(npsiList)==1 && length(npsiList[[1]])== length(nomiPS) ) {
      id.ps <- seq(id.ps, l=length(nomiPS))
      npsiList<-sapply(npsiList[[1]], as.list)
      estList <- rep(estList, length(nomiPS))
      psiList <- rep(psiList, length(nomiPS))
      byList <- rep(byList, length(nomiPS) )
      nomiBy <-  rep(nomiBy, length(nomiPS) )
      RList <- rep(RList, length(nomiPS) )
      fixpsiList <- rep(fixpsiList, length(nomiPS) )
      fxList <- rep(fxList, length(nomiPS) )
      levelsBy <-rep(levelsBy, length(nomiPS) )
      mf<- as.data.frame(sapply(1:ncol(mf), function(.x) mf[.x]))
      names(mf)[id.ps]<-nomiPS
    }
    
    
    
    if(length(nomiPS)!=length(npsiList)){ #prima era length(psiList), ma dovrebbe essere lo stesso.. L'ho cambiato per tener conto di X matrice..
      id.ps <- seq(id.ps, l=length(nomiPS))
      npsiList <- rep(npsiList, length(nomiPS) )
      estList <- rep(estList, length(nomiPS))
      psiList <- rep(psiList, length(nomiPS))
      byList <- rep(byList, length(nomiPS) )
      nomiBy <-  rep(nomiBy, length(nomiPS) )
      RList <- rep(RList, length(nomiPS) )
      fixpsiList <- rep(fixpsiList, length(nomiPS) )
      fxList <- rep(fxList, length(nomiPS) )
      levelsBy <-rep(levelsBy, length(nomiPS) )
      mf<- as.data.frame(sapply(1:ncol(mf), function(.x) mf[.x]))
      names(mf)[id.ps]<-nomiPS
    }
    
    if(all(sapply(levelsBy, is.null)) && (length(npsiList)!=length(nomiPS))) stop(" 'npsi' is not correctly specified")
    
    
    #testo.ps <-names(mf)[id.ps]
    nomiCoefUNPEN<-names(mf)[-c(1,id.ps)]
    X <- if(!is.empty.model(mt)){
      model.matrix(mt, mf, contrasts) 
    } else {stop("error in the design matrix")}#matrix(, NROW(Y), 0L)
    attrContr<-attr(X, "contrasts")
    n<-nrow(X)
    
    
    
    limZ<- rangeSmooth<-mVariabili<-NULL
		#se ci sono termini ps()+ps(..,by) il nome delle variabili smooth vengono cambiati per aggiungere la variabile by
		nomiPS.orig <- nomiPS
		  
		nomiPS.By <-paste(nomiPS, nomiBy, sep=":")
		nomiPS <- unlist(lapply(nomiPS.By, function(.x) sub(":NULL", "", .x)))

		####se la stessa variabile e' specificata come ps o termine lineare..
		#if(length(intersect(nomiPS,nomiCoefUNPEN))>=1) stop("The same variable specified in seg() and as linear term")

		#ATTENZIONE.. se vuoi subito costruire i nomi ps(x), ps(x):z, ecc...usa i:
		nomiPS.ps<- sapply(nomiPS.orig, function(.x)paste("seg(",list(.x),")",sep=""))
		nomiPS.ps<-unlist(lapply(paste(nomiPS.ps, nomiBy, sep=":"), function(.x) sub(":NULL", "", .x)))
		nomiPS.ps.list<-as.list(nomiPS.ps) #serve lista


		for(j in id.ps) mVariabili[length(mVariabili)+1]<-mf[j]
    B<- Bfix <- nomiPS.ps.int.list<- nomiPS.ps.int.list.All<- vector(length=length(mVariabili) , "list")
    #BFixed<-BB<-Bderiv
    
    #browser()
    variabiliBY<-vector("list", length(mVariabili))
    for(j in 1:length(mVariabili)) {
      if(nomiBy[j]=="NULL"){ # se usuale termine seg()
        nomiBy.j <- NULL
        variabileSmooth<- c(mVariabili[[j]]) #c() converte le matrici in vettori, drop() no..!
        variabileSmooth<- fxList[[j]](variabileSmooth)
        #for(jj in c("nomeX", "psi", "npsi", "f.x", "nomeBy")) attr(variabileSmooth,jj)<-NULL
        B[[j]]<- variabileSmooth
        rangeSmooth[[j]] <- range(variabileSmooth)
        limZ[[j]] <- quantile(variabileSmooth, names=FALSE, probs=c(alpha[1],alpha[2]))
        nomiPS.ps.int.list[[j]]<- nomiPS[j]
        nomiPS.ps.int.list.All[[j]]<- nomiPS.orig[j]
        
      } else { #se ci sono termini by
        if(is.null(levelsBy[[j]])){ #se e' vc con variabile continua
          
          nomiBy.j <- nomiBy[j]
          B[[j]] <- variabileSmooth<- mVariabili[[j]][,1,drop=TRUE] 
          #variabiliBY[[j]] <- mVariabili[[j]][, -1,drop=TRUE]
          #byList include le variabili by..
          #nomiCoefPEN[[j]]<- sapply(1:ncol(B[[j]]), function(x) gsub(":", paste(".",x, ":", sep="") , nomiPS.ps[j]))
          rangeSmooth[[j]] <- range(variabileSmooth)
          limZ[[j]] <- quantile(variabileSmooth, names=FALSE, probs=c(alpha[1],alpha[2]))
          nomiPS.ps.int.list[[j]]<- nomiPS[j]
          nomiPS.ps.int.list.All[[j]]<- paste(nomiPS.orig[j], strsplit(nomiBy[j], ",")[[1]], sep=":")
          if(any(strsplit(nomiBy[j], ",")[[1]]=="")){ #se c'e' un termine senza nome (seg(tt, by=cbind(1,x)))
            .id <- which(strsplit(nomiBy[j], ",")[[1]]=="")
            for(.idj in 1:length(.id)) nomiPS.ps.int.list.All[[j]][.idj] <- gsub(":", "", nomiPS.ps.int.list.All[[j]][.idj])
          }
        } else {#se e' VC con variabile categoriale
          nomiBy.j <- nomiBy[j]
          variabileSmooth<- mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE] 
          variabileSmooth<- fxList[[j]](variabileSmooth) #attr( mVariabili[[j]], "f.x")(variabileSmooth)
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
          nomiPS.ps.int.list.All[[j]]<- nomiPS.ps.int.list[[j]]
          
        }
      } 
    } #end for(j in 1:length(mVariabili))
    allNameOK <- unlist(nomiPS.ps.int.list.All)
    if(length(allNameOK)!=length(unique(allNameOK))) stop("some error in specification of seg() terms")
    
    #browser()
    
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
      
    if(!any(sapply(psiList,is.list))) psiList <- rep(psiList, repl)
    if(!any(sapply(estList,is.list))) estList <- rep(estList, repl)
    if(!any(sapply(RList,is.list))) RList <- rep(RList, repl)
    nomiPS.orig <- rep(nomiPS.orig, repl)
    Bfix <- rep(Bfix, repl)
    fixpsiList<- rep(fixpsiList, repl)
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
        
        #se la lista contiene solo NULL, non funziona...
        #penMatrixList <- append(penMatrixList, penMatrixList[[id.vc]], after = id.vc-1)
        #penMatrixList[[id.vc+nc]]<-NULL
    }
    
    #browser()
    
    
    
      #if(!all(sapply(estList,check.estPsi))) stop(" 'est' is misspecified in one or more seg() term")
      
    nomiTerminiSEG<-nomiCoefPSI <- nomiCoefU <- NULL
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
    for(j in 1:length(B)){
        K <- if(!is.na(npsiList[nomiSeg[j]])) npsiList[nomiSeg[j]] else npsiList[j]
        npsiList1[j]<- K
        if(any(is.na(psiList[[j]]))){ 
          if(control$quant) {
            psiList[[j]]<- quantile(B[[j]], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)
          } else {
            psiList[[j]]<- (min(B[[j]])+ diff(range(B[[j]]))*(1:K)/(K+1))
          }
        } else {
          K<-npsiList1[j]<-length(psiList[[j]])
        }
        if(!is.null(fixpsiList[[j]])) {
          Bfix[[j]]<- sapply(sort(fixpsiList[[j]]), function(.x) 1*(B[[j]]>.x))
          colnames(Bfix[[j]])<- paste("U", 1:length(fixpsiList[[j]]),".fixed.",nomiPS.orig[j], sep="")
        }
        #se per qualche termine ci sono le matrici dei vincoli sulle slope
        j.ok=match(nomiSeg[j], names(RList), nomatch=0)
        j.ok <-if(j.ok>0) j.ok else j 
        if(!any(is.na(RList[[j.ok]]))){
          RList[[j]] <- RList[[j.ok]]
          id.contrR[j] <-TRUE
          } else {
            j.ok=match(nomiSeg[j], names(estList), nomatch=0)
            j.ok <-if(j.ok>0) j.ok else j
            if(!any(is.na(estList[[j.ok]]))){
              if(length(estList[[j.ok]])!=(K+1)) stop(" 'est' is not compatible with 'n.psi' ")
              #browser()
              RList[[j]]<-diag(K+1)[,estList[[j.ok]]==1,drop=FALSE]
              id.contrR[j] <-TRUE
            } else {
              RList[[j]]<-diag(K+1)
              id.contrR[j] <-FALSE
            }
          }
        #nomiCoefPSI[[j]]<- paste(paste("psi",1:length(psiList[[j]]), sep=""), nomiPS.ps.unlist[j], sep=".") 
        #nomiCoefU[[j]]<- paste(paste("U",1:length(psiList[[j]]), sep=""), nomiPS.ps.int.list.All[[j]], sep=".")
        #nomiTerminiSEG[[j]]<-rep(nomiPS.ps.unlist[j], length(psiList[[j]]))
      }
    
    for(i in 1:length(B)) {
      nomiCoefPSI[[i]]<- paste(paste("psi",1:length(psiList[[i]]), sep=""), nomiPS.ps.unlist[i], sep=".") ##oppure sep=".psi"
      nomiTerminiSEG[[i]]<-rep(nomiPS.ps.unlist[i], length(psiList[[i]]))
    }
    #browser()
    
    nomiCoefU<-lapply(nomiCoefPSI, function(.x) sub("psi","U",.x )) 
    nomiCoefZ<-lapply(nomiCoefPSI, function(.x) sub("psi","Z",.x ))
    
    #

      npsii <- sapply(psiList,length)
      id.psi.group <- rep(1:length(psiList), npsii)
      Z<- lapply(1:length(B), function(.x) matrix(B[[.x]], nrow=n, ncol=npsiList1[[.x]]))
      Z<- do.call(cbind,Z)
      
      #NB11111 colnames(Z) <- unlist(nomiCoefZ) #unlist(nomiTerminiSEG)
      colnames(Z) <- unlist(nomiTerminiSEG)
      
      
      
      
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
      #browser()
      X<- X[, !startsWith(colnames(X),"seg("), drop=FALSE]

      idZ <- unlist(tapply(id.psi.group, id.psi.group, function(.x) c(TRUE, rep(FALSE, length(.x)-1))))
      Z.ok<-Z[, idZ, drop=FALSE]
      colnames(Z.ok) <- nomiPS.ps.unlist
      
      #X<-cbind(X, Z.ok) #Z.ok include anche i termini lineari delle variabili segmented
      
      #colnames(Z)<- unlist(nomiCoefPEN)
      initial <- unlist(psiList)
      PSI <- matrix(initial, n, length(initial), byrow = TRUE)
      
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
      if(is.null(weights)) weights<-rep(1,n)
      orig.offs<-offs
      if(is.null(offs)) offs<-rep(0,n)

      invXtX<-Xty<-NULL
      if(is.null(alpha)) alpha<- max(.05, 1/length(Y))
      if(length(alpha)==1) alpha<-c(alpha, 1-alpha)
      
      #rangeZ <- apply(Z, 2, range) #Z ha molti 0 se e' il prodotto con dummy (quando c'e' seg(x, by))
      
      limZ<-do.call(cbind, lapply(limZ, function(.x){if(is.list(.x)) do.call(cbind, .x) else cbind(.x)} ))
      #limZ<-matrix(sapply(1:length(npsii), function(.x) rep(limZ[,.x],npsii[.x])), nrow=2, byrow = FALSE)
      limZ <- do.call(cbind, lapply(1:length(npsii), function(.x) matrix(limZ[,.x],nrow=2,ncol=npsii[.x])))
      rangeZ <- do.call(cbind, lapply(1:length(npsii), function(.x) matrix(rangeSmooth[[.x]],nrow=2,ncol=npsii[.x])))
      #browser()
      colnames(rangeZ) <- unlist(nomiTerminiSEG)

      
      #sapply(byList, function(.x) {if(length(.x)>0 && !is.matrix(.x)) 1 else ncol(.x)})
      #browser()
      #13/03/24: ho tolto dev0=var(Y)*(n-1) perche'
      #   con Pois con contegg bassi viene molto piccola!! e quindi l'algoritmo non partiva!
      #   e comunque si dovrebbe chiamare L0
      #
      #browser()
      byList <- lapply(byList, function(.x) if(is.vector(.x)) matrix(.x) else .x)
      opz<-list(toll=toll,h=h,stop.if.error=stop.if.error, display=display,it.max=it.max,nomiOK=unlist(nomiCoefU), usestepreg=TRUE,
                fam=family, eta0=NULL, maxit.glm=maxit.glm, id.psi.group=id.psi.group, gap=gap, limZ=limZ,
                conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step,
                pow=pow, #visualBoot=visualBoot, 
                digits=digits, fc=fc, RList=RList, nomiSeg=nomiSeg, seed=control$seed,
                npsii=npsii, agg=agg, byList=byList,rangeZ=rangeZ, tol.opt=control$tol.opt)
      
      #browser()
      
      if(any(sapply(levelsBy, is.null)) && any(!sapply(byList, is.null))){ #se ci sono Struct Changes
        idSC<-TRUE
        
        if(fitter0=="lm"){
        if(n.boot <=0 ) {
          obj <- step.lm.fitSC(Y, X, Z, PSI, weights, offs, opz)
          return(obj)
        } else {
          obj <- step.lm.fitSC.boot(Y, X, Z, PSI, weights, offs, opz, 
                                  n.boot = n.boot, size.boot = size.boot, random = random, 
                                  break.boot = break.boot)
          seed<- obj$seed
        }
          } else {
            if(n.boot <=0 ) {
              obj <- step.glm.fitSC(Y, X, Z, PSI, weights, offs, opz)
              return(obj)
              } else {
                obj <- step.glm.fitSC.boot(Y, X, Z, PSI, weights, offs, opz, 
                                      n.boot = n.boot, size.boot = size.boot, random = random, 
                                      break.boot = break.boot)
                seed<- obj$seed
              }
          }
        } else {
          idSC<-FALSE
          if(any(id.contrR)){
            if(fitter0=="lm"){
              if(n.boot <= 0) {
                obj <- stepConstr.lm.fit(Y, X, Z, PSI, weights, offs, opz)
          } else {
            obj <- stepConstr.lm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                                 n.boot = n.boot, size.boot = size.boot, random = random, 
                                 break.boot = break.boot)
            seed<- obj$seed
          }
          class0<- "lm"
          if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
        } else {
          if(n.boot<=0){
            obj <-stepConstr.glm.fit(Y, X, Z, PSI, weights, offs, opz)
          } else {
            obj <-stepConstr.glm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                                   n.boot=n.boot, size.boot=size.boot, random=random, 
                                   break.boot=break.boot)
            seed<- obj$seed
          }
          class0<-c("glm","lm")
          }
      } else {
        if(fitter0=="lm"){
          if(n.boot <= 0) {
            obj <- step.lm.fit(Y, X, Z, PSI, weights, offs, opz)
            } else {
            #browser()
              obj <- step.lm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                               n.boot = n.boot, size.boot = size.boot, random = random, 
                               break.boot = break.boot)
              seed<- obj$seed
          }
          class0<-"lm"
          #if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
        } else {
          if(n.boot<=0){
            obj <-step.glm.fit(Y, X, Z, PSI, weights, offs, opz)
          } else {
            obj <-step.glm.fit.boot(Y, X, Z, PSI, weights, offs, opz, 
                                n.boot=n.boot, size.boot=size.boot, random=random, 
                                break.boot=break.boot)
            seed<- obj$seed
        }
        class0<-c("glm","lm")
        }
      }
      }
      
      
      #browser()
      
      #da modificare... vedi stepmented.lm o stepmented.glm
      
      if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(invisible(NULL))
      }
      #
      id.warn <- obj$id.warn
      it <- obj$it
      epsilon <- obj$epsilon
      psi<-obj$psi
      psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
      #id.psi.group<-obj$id.psi.group #id.psi.group c'e' gia'....
      npsi.groups <- tapply(id.psi.group, id.psi.group, length)
      
      #Nelle funzioni step servono solo i psi e le U (i beta.c NON servono)
      #i psi sono gia' ordinati
      #psi<-unlist(tapply(psi, id.psi.group, sort)) 
      Z0 <-apply(Z,2,sort)
      npsi<- sum(npsii)
      
      ris.psi<-cbind(Est.=psi, St.Err=NA)
      name.Z <- unlist(nomiTerminiSEG) #e nel caso di Structural changes?
      U <- obj$U #le U sono calcolate sui psi non-rounded, pero' mi sa che e' lo stesso...
      
      psi.rounded<-sapply(1:npsi, function(j) Z0[sum(Z0[,j]<psi[j])+c(0,1,2),j])
      #psi.rounded<-sapply(1:npsi, function(j) Z0[sum(Z0[,j]<psi[j])+c(-1,0,1,2),j]) #se vuoi prendere il precedente..
      #browser()

      id.new.result<-rep(FALSE, npsi)
      
      if(control$check.next){ #se vuoi fare un controllo sulle soluzioni +1..
        L0 <- as.numeric(obj$SumSquares.no.gap)
        ############== definisci la f
        final.fit.f<-function(fitter0){
          if(fitter0=="lm"){
            if(is.null(weights)){
              objV <- .lm.fit(cbind(X, U), (Y- offs))
              r<- sum(objV$residuals^2)
            } else {
              sw <- sqrt(weights)
              objV <- .lm.fit(sw*cbind(X, U), sw*(Y- offs))
              r<-sum(weights*objV$residuals^2)
              }
            } else {#if glm
              eta0 <- attr(obj$SumSquares.no.gap, "eta") #obj$eta0
              objV <- try(suppressWarnings(glm.fit(cbind(X, U), y = Y, offset = offs,
                                                 weights = weights, family = opz$fam, #control = glm.control(maxit = maxit.glm), 
                                                 etastart = eta0)), silent = TRUE) #obj$obj$linear.predictors
              r<-objV$deviance
            }
          r
        }
        psi.roundedOK <- psi.rounded[1:2,,drop=FALSE]
        psiTry<-psi.rounded[1,]
        
        if(idSC) {
          PSI<- matrix(psi.rounded[1,], n, npsi, byrow = TRUE)
          idUpsi <- rep(1:npsi, sapply(byList, function(.x) if(is.null(.x)) 1 else ncol(.x )))#per individuare le U corrispondenti ai 
          #diversi psi.. Infatti con SC si possono avere diverse U per un unico psi
          for(j in 1:npsi){
            psiTry[j] <- psi.rounded[2,j]
            Unew <- obj$fn.U(Z, matrix(psiTry, n, npsi, byrow = TRUE), byList, id.psi.group)
            U[,idUpsi==j] <- Unew[,idUpsi==j]
            Lnew <- final.fit.f(fitter0)
            if(Lnew<L0) {#se il fit e' migliore prendi quella soluzione... e lascia la U che e' gia' calcolata per quella soluzione
              psi.roundedOK[,j]<- psi.rounded[-1,j]
              id.new.result[j]<-TRUE
            } else { #.. altrimenti ri-calcola la U con la soluzione precedente..
                Unew <- obj$fn.U(Z, PSI, byList, id.psi.group)
                U[,idUpsi==j] <- Unew[,idUpsi==j]
            }
          }
        } else {
          #browser()
          for(j in 1:npsi){
            psiTry[j] <- psi.rounded[2,j]
            U[,j]<- 1*(Z[,j]>psiTry[j])
            Lnew <- final.fit.f(fitter0)
            if(Lnew<L0) {#se il fit e' migliore prendi quella soluzione... e lascia la U che e' gia' calcolata per quella soluzione
                psi.roundedOK[,j]<- psi.rounded[-1,j]
                id.new.result[j]<-TRUE
              } else { #.. altrimenti ri-calcola la U con la soluzione precedente..
                U[,j]<- 1*(Z[,j]>psi.rounded[1,j])
              }
            }
        }
        psi.rounded <- psi.roundedOK
        
      }
        
      
      #se vuoi fare un controllo anche sulla soluzione precedente
      #questo potrebbe servire per aumentare un po' la toll di optimize()???
      
      # psi.rounded<-sapply(1:npsi, function(j) Z0[sum(Z0[,j]<psi[j])+c(-1,0,1,2),j])
      # L0 <- as.numeric(obj$SumSquares.no.gap)
      # psi.roundedOK <- psi.rounded[2:3,,drop=FALSE]
      # 
      # psiTry<-psiTry1<-psiTry2<-psi.rounded[2,]
      # 
      # for(j in 1:npsi){
      #   psiTry1[j] <- psi.rounded[1,j] #il precedente
      #   U[,j]<- 1*(Z[,j]>psiTry1[j])
      #   Lnew1 <- final.fit.f(fitter0)
      #   psiTry2[j] <- psi.rounded[3,j] #il successivo
      #   U[,j]<- 1*(Z[,j]>psiTry2[j])
      #   Lnew2 <- final.fit.f(fitter0)
      #   
      #   psiOK.j <- c(psiTry1[j], psi.roundedOK[2,j] , psiTry2[j])[which.min(Lnew1, L0, Lnew2)]
      #   U[,j]<- 1*(Z[,j]>psiOK.j)
      #   psi.roundedOK[1,j]<-psiOK.j
      # }
      # 
      # all.psi.ok<-psi.roundedOK[1,] #include le soluzioni
      # #il problema e' 
      # id=sapply(1:ncol(psi.rounded), function(.x) which(all.psi.ok[.x]==data.frame(psi.rounded)[[.x]]))
      # psi.roundedOK<-sapply(1:ncol(psi.rounded), function(.x) psi.rounded[id[.x]+0:1,.x])
      # 


      colnames(U)<- nomiU <-unlist(nomiCoefU)
      nomiVxb <- unlist(nomiCoefPSI)
      nomiV<- gsub("psi", "V", nomiVxb)
      colnames(psi.rounded)<-rownames(ris.psi)<-names(psi)<-nomiVxb
      rownames(psi.rounded)<-c("inf [","sup (")

      se.psi<-rep(NA, npsi)
      if(fitter0=="lm"){
        class0<-"lm"
        objV <- if(is.null(weights)) lm.fit(cbind(X, U), Y, offset = offs) else lm.wfit(cbind(X, U), Y, weights, offset = offs)
        objV$df.residual <- objV$df.residual- length(psi)
        L0 <- sum(weights*objV$residuals^2)
        if(var.psi) {
          s2 <- L0/objV$df.residual
          R <- chol2inv(objV$qr$qr)
          se.psi <- sqrt(diag(R)*s2)[match(nomiVxb, names(coef(objV)),0)]
        }
      } else {
        class0<-c("glm", "lm")
        eta0 <- attr(obj$SumSquares.no.gap, "eta") #obj$eta0
        objV <- try(suppressWarnings(glm.fit(cbind(X, U), y = Y, offset = offs,
                                             weights = weights, family = opz$fam, #control = glm.control(maxit = maxit.glm), 
                                             etastart = eta0)), silent = TRUE) #obj$obj$linear.predictors
        objV$df.residual <- objV$df.residual- length(psi)
        L0 <- objV$deviance
        if (length(offs) && attr(mt, "intercept") > 0L) {
          #se c'e' un offset devi calcolare la null.deviance (come fa glm())
          obj0 <- try(suppressWarnings(glm.fit(X[, "(Intercept)", drop = FALSE], y = Y, offset = offs,
                                               weights = weights, family = opz$fam, #control = glm.control(maxit = maxit.glm), 
                                               etastart = eta0, intercept=TRUE)), silent = TRUE)
          
          # obj0 <- eval(call(if (is.function(method)) "method" else method, 
          #                   x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values, 
          #                   weights = weights, offset = offset, family = family, 
          #                   control = control, intercept = TRUE))
          if (!obj0$converged) warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
          objV$null.deviance <- obj0$deviance
        }
        
        if(var.psi) {
          R <- chol2inv(objV$qr$qr)
          s2 <- 1
          if(!opz$fam$fam%in%c("poisson","binomial")) s2<- L0/objV$df.residual
          se.psi <- sqrt(diag(R)*s2)[match(nomiVxb, names(coef(objV)),0)]
        }
      }
      
      
      if(any(id.new.result) && display) cat(" Better objective found:", L0, "at psi =", psi.rounded[1,], "\n")
      
      objV$rank <- objV$rank + length(psi)
      if(!is.null(objV$aic)){
        objV$aic <- objV$aic +  2*length(psi)
      }
      #objV$nameUV <- list(U = drop(nomiU), V = nomiV, Z = name.Z) #Z = name.Z
      #browser()
      #objV$nameUV$formulaSegAllTerms<- paste("~",paste(nomiSeg, collapse="+"))
      #paste("~", paste(sapply(all.seg.form, function(.x) strsplit(paste(.x), "~"))[2,],collapse="+"))
      
      objV$rangeZ<-obj$rangeZ
      objV$call <- match.call()
      objV$psi<-ris.psi
      objV$psi.history <- psi.values
      objV$psi.rounded <- psi.rounded
      if(n.boot>0) objV$seed <- seed
      
      
      #browser()
      
      objV$Z <- Z.ok #Z[,unique(name.Z),drop=FALSE]
      
      
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
      
      names(all.seg.form)<-nomiPS
      #costruisci la formulaLin.. Attenzione non tiene conto di eventuali vincoli sulle pendenze.
      splitFo <- strsplit(as.character(formula),"[+]")
      #allX.lin<-paste(c(splitFo[[3]][-grep("seg[(]", splitFo[[3]])], unique(nomiPS.orig)), collapse="+") #anche le variabili seg
      
      termLin <- splitFo[[3]][-grep("seg[(]", splitFo[[3]])]
      if(length(termLin)>0){ #se ci sono altre variabili lineari (i.e. non-seg)
        allX.lin <- paste(termLin , collapse="+")                          #solo i termini non-seg!
        formulaLin <- as.formula(paste(splitFo[[2]], splitFo[[1]], allX.lin)) #formula escluso i termini seg 
        nomiVarLin<- setdiff(all.vars(formulaLin), splitFo[[2]]) #vettore di nomi ("x", "z" )
        #termLin <- strsplit(allX.lin,"[+]")[[1]] #se allX.lin e' una formula, vettore di termini ("poly(x,2)", "z")
        Z.in.obj<-intersect(nomiVarLin, nomiSeg) #nomiSeg=name.Z
        
        if(length(Z.in.obj)>0){
          f.x<-matrix(NA, 150, ncol(Z.ok[,Z.in.obj,drop=FALSE])) #prima era nrow(objF$Z) invece che 100
          for(j in 1:length(Z.in.obj)){
            termLin.ok <- grep(Z.in.obj[j], termLin, value=TRUE)
            dd<-data.frame(seq(min(Z.ok[,Z.in.obj[j]]), max(Z.ok[,Z.in.obj[j]]), l=nrow(f.x)))
            names(dd)<- Z.in.obj[j]
            M <- model.matrix(reformulate(termLin.ok, intercept=FALSE), data=dd)
            f.x[,j]<-M%*% coef(objV)[colnames(M)]
          }
          colnames(f.x)<-Z.in.obj
          objV$f.x<-f.x
        }
      } else {
        formulaLin <- if(interc) as.formula(paste(splitFo[[2]], splitFo[[1]], "1")) else as.formula(paste(splitFo[[2]], splitFo[[1]], "0")) 
      }
      
      #browser()
      
      objV$nameUV <- list(U = drop(nomiU), V = nomiV, Z = name.Z)
      #objV$nameUV <- list(U = drop(nomiU), V = rownames(ris.psi), Z = nomiPS.orig) #nomiPS.orig??
      objV$nameUV$formulaSeg<- all.seg.form
      objV$nameUV$formulaSegAllTerms<- paste("~", paste(sapply(all.seg.form, function(.x) strsplit(paste(.x), "~"))[2,],collapse="+"))
      
      objV$formulaLin<- formulaLin
      objV$terms <- mt
      objV$y<-Y
      if(x) objV$x <- X
      objV$contrasts <- attrContr  
      objV$xlevels <- .xlivelli
      objV$it <- it 
      objV$epsilon <- epsilon
      objV$id.warn <- id.warn
      objV<- structure(c(objV, list(offset=orig.offs)))
      class(objV)<-c("stepmented", class0)
      objV
}
      
      
      
     
      
      
