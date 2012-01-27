intercept<-function (ogg, parm, gap=TRUE, rev.sgn = FALSE, var.diff = FALSE){
#E' come sotto ma unisce due cicli for..
#corregge in caso di no model intercept
#gap: mostra le interc con gap ..
    if (!"segmented" %in% class(ogg))
        stop("A segmented model is needed")
    if (var.diff && length(ogg$nameUV$Z) > 1) {
        var.diff <- FALSE
        warning("var.diff set to FALSE with multiple segmented variables",
            call. = FALSE)
    }
    nomepsi <- rownames(ogg$psi)
    nomeU <- ogg$nameUV[[1]]
    nomeZ <- ogg$nameUV[[3]]
    if (missing(parm)) {
        nomeZ <- ogg$nameUV[[3]]
        if (length(rev.sgn) == 1)
            rev.sgn <- rep(rev.sgn, length(nomeZ))
    }
    else {
        if (!all(parm %in% ogg$nameUV[[3]])) {
            stop("invalid parm")
        }
        else {
            nomeZ <- parm
        }
    }
    if (length(rev.sgn) != length(nomeZ))
        rev.sgn <- rep(rev.sgn, length.out = length(nomeZ))
    nomi <- names(coef(ogg))
    nomi <- nomi[-match(nomepsi, nomi)]
    Allpsi <- index <- vector(mode = "list", length = length(nomeZ))
    gapCoef<-summary.segmented(ogg)$gap
    Ris <- list()
    digits <- max(3, getOption("digits") - 3)
    rev.sgn <- rep(rev.sgn, length.out = length(nomeZ))
    if("(Intercept)"%in%names(coef(ogg))){
        alpha0 <- alpha00 <- coef(ogg)["(Intercept)"]} else {alpha0 <- alpha00 <-0}
#per ogni variabile segmented...
    for (i in 1:length(nomeZ)) {
        id.cof.U <- grep(paste("\\.", nomeZ[i], "$", sep = ""),
            nomi, value = FALSE)
        psii <- ogg$psi[grep(paste("\\.", nomeZ[i], "$", sep = ""),
            rownames(ogg$psi), value = FALSE), 2]
        Allpsi[[i]] <- psii
        id.cof.U <- id.cof.U[order(psii)]
        index[[i]] <- id.cof.U
        alpha0<-if("(Intercept)"%in%names(coef(ogg))) coef(ogg)["(Intercept)"] else 0
        ind <- as.numeric(na.omit(unlist(index[[i]])))
        cof <- coef(ogg)[ind]
        alpha <- vector(length = length(ind))
        gapCoef.i<-gapCoef[grep(paste("\\.",nomeZ[i],"$",sep=""), rownames(gapCoef), value = FALSE),"Est."]
        for (j in 1:length(cof)) {
            alpha[j] <- alpha0 - Allpsi[[i]][j] * cof[j]
            if(gap) alpha[j] <- alpha[j] - gapCoef.i[j]
            alpha0 <- alpha[j]

        }
        #if(gap) alpha<-alpha -gapCoef[grep(paste("\\.",nomeZ[i],"$",sep=""), rownames(gapCoef), value = FALSE),"Est."]
        cof.out <- c(alpha00, alpha)
        ris <- matrix(cof.out)
        dimnames(ris) <- list(paste("intercept", 1:nrow(ris),
            sep = ""), "Est.")
        Ris[[nomeZ[i]]] <- signif(ris, digits)
    }
    Ris
}
