intercept<-function (ogg, parm, rev.sgn = FALSE, var.diff = FALSE, .vcov=NULL, .coef=NULL,
    digits = max(4, getOption("digits") - 2),...){
#corregge in caso di no model intercept -- CHE VOLEVO DIRE?? #forse che adesso funziona se nel modello non c'e' l'interc.
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
    
    covv <- if(is.null(.vcov)) vcov(ogg, ...) else .vcov 
    estcoef<- if(is.null(.coef)) coef(ogg) else .coef
    
    if(!all(dim(covv)==c(length(estcoef), length(estcoef)))) stop("dimension of cov matrix and estimated coeffs do not match", call. = FALSE)
    
    
    if (var.diff && length(ogg$nameUV$Z) > 1) {
        var.diff <- FALSE
        warning("var.diff set to FALSE with multiple segmented variables",
            call. = FALSE)
    }
    #browser()
    nomepsi <- rownames(ogg$psi)
    nomeU <- ogg$nameUV$U
    nomeZ <- ogg$nameUV$Z
    if (missing(parm)) {
        nomeZ <- ogg$nameUV$Z
        if (length(rev.sgn) == 1)
            rev.sgn <- rep(rev.sgn, length(nomeZ))
    } else {
        if (!all(parm %in%nomeZ)) {
            stop("invalid parm")
        } else {
            nomeZ <- parm
        }
    }
    if (length(rev.sgn) != length(nomeZ)) rev.sgn <- rep(rev.sgn, length.out = length(nomeZ))
    nomi <- names(estcoef)
    #browser()
    nomi <- nomi[-match(nomepsi, nomi)]
    Allpsi <- index <- vector(mode = "list", length = length(nomeZ))
    # gapCoef<-summary.segmented(ogg)$gap   ##eliminato 10/11/15
    Ris <- list()
    rev.sgn <- rep(rev.sgn, length.out = length(nomeZ))
    
    alpha0 <- alpha00 <-0
    idInterc<-grep("ntercept",names(estcoef))
    if(length(idInterc)>0) if(idInterc!= grep("ntercept",rownames(covv))) stop("intercept name in coeff and vcov do not match")
    if(length(idInterc)==1) alpha0 <- alpha00 <- estcoef[idInterc]
    
#    if("(Intercept)"%in%names(estcoef) || "intercept"%in%names(estcoef)){
#        alpha0 <- alpha00 <- estcoef["(Intercept)"]} else {alpha0 <- alpha00 <-0}
    #per ogni variabile segmented...
    #browser()
    for (i in 1:length(nomeZ)) {
        #id.cof.U <- f.U(ogg$nameUV$U, nomeZ[i]) + (match(ogg$nameUV$U[1], nomi)-1)
        #psii<- ogg$psi[f.U(ogg$nameUV$V, nomeZ[i]) , "Est."]
        #Allpsi[[i]] <- sort(psii, decreasing = FALSE) 
        #id.cof.U <- id.cof.U[order(psii)]
        #index[[i]] <- id.cof.U
        #ind <- as.numeric(na.omit(unlist(index[[i]])))
        #cof <- coef(ogg)[ind]
        #alpha0<-if("(Intercept)"%in%names(estcoef)) estcoef["(Intercept)"] else 0
        alpha0<-if(length(idInterc)==1) estcoef[idInterc] else 0

        Allpsi[[i]] <-ogg$indexU[[nomeZ[i]]]
        
        if(is.null(ogg$constr)){
          cof<- estcoef[names(ogg$indexU[[nomeZ[i]]])]
        } else {
          index <- match(c(nomeZ[i],ogg$nameUV$U[grep(nomeZ[i], ogg$nameUV$U)]), names(coef(ogg)),0)
          cof<- drop(ogg$constr$invA.RList[[match(nomeZ[i], ogg$nameUV$Z,0)]]%*%coef(ogg)[index])[-1] #solo le 
        }
        alpha <- vector(length = length(cof)) #length(ind)
        for (j in 1:length(cof)) {
            alpha[j] <- alpha0 - Allpsi[[i]][j] * cof[j]
            alpha0 <- alpha[j]
        }
        cof.out <- c(alpha00, alpha)
        if(rev.sgn[i]) cof.out <- cof.out[length(cof.out):1]
        ris <- matrix(cof.out)
        dimnames(ris) <- list(paste("intercept", 1:nrow(ris), sep = ""), "Est.")
        Ris[[nomeZ[i]]] <- signif(ris, digits)
    }
    Ris
}
