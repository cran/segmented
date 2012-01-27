broken.line<-function (ogg, term = NULL, gap = FALSE, linkinv = FALSE, interc = TRUE)
{
    if (!"segmented" %in% class(ogg))
        stop("A segmented model is requested")
    nomepsi <- rownames(ogg$psi)
    nomeU <- ogg$nameUV$U
    nomeZ <- ogg$nameUV$Z
    nomiSenzaV <- nomiSenzaU <- nomi <- names(coef(ogg))
    nomiSenzaU[match(nomeU, nomi)] <- ""
    nomiSenzaV[match(nomepsi, nomi)] <- ""
    index <- vector(mode = "list", length = length(nomeZ))
    for (i in 1:length(nomeZ)) {
        index[[i]] <- c(match(nomeZ[i], nomi), grep(paste("\\.",
            nomeZ[i], "$", sep = ""), nomiSenzaV, value = FALSE))
        if (gap)
            index[[i]] <- c(index[[i]], grep(paste("\\.", nomeZ[i],
                "$", sep = ""), nomiSenzaU, value = FALSE))
    }
    variabili <- Ris <- list()
    for (i in 1:length(index)) {
        ind <- as.numeric(na.omit(unlist(index[[i]])))
        cof <- coef(ogg)[ind]
        Ris[[nomeZ[i]]] <- cof
        variabili[[nomeZ[i]]] <- data.matrix(ogg$model[names(cof)])
    }
    ris <- mapply(function(xx, yy) drop(xx %*% yy), variabili,Ris)
    if (interc)
        ris <- ris + coef(ogg)["(Intercept)"]
    if (!is.null(term))
        ris <- ris[, term, drop=FALSE]
    if (inherits(ogg, what = "glm", FALSE) && linkinv)
        ris <- apply(ris, 2, ogg$family$linkinv)
    return(ris)
}
