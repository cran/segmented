`davies.test` <- 
function (ogg, term, k = 10, alternative = c("two.sided", "less", "greater")) {
#ho messo as.expression(formula(ogg)) potresti anche as.expression(ogg$call)
#Non funziona se ogg non è chiamato con l'argomento data.
    if (!is.character(term))
        stop("term must be character")
    alternative <- match.arg(alternative)
    if(!term%in%all.names(formula(ogg))) stop("term has be included in the model")
    Z <- if (term %in% names(ogg$model)) {
        ogg$model[, term]
    }
    else {
        eval(parse(text = term))
    }
    qq <- quantile(Z, prob = c(0.05, 0.95), names = FALSE, na.rm = TRUE)
    sx <- qq[1]
    dx <- qq[2]
    valori <- seq(sx, dx, length = k)
    ris.valori <- NULL
    for (i in valori) {
        o <- segmented(ogg, seg.Z = ~Z, psi = list(Z = i), control = seg.control(it.max = 0,
            display = FALSE))
        ris.valori[(length(ris.valori)) + 1] <- if (is.list(o)) {
            summary(o)$coef[length(coef(ogg)) + 1, 3]
        }
        else {
            NA
        }
    }
    valori <- valori[!is.na(ris.valori)]
    ris.valori <- ris.valori[!is.na(ris.valori)]
    V <- sum(abs(diff(ris.valori)))
    onesided <- TRUE
    if (alternative == "less") {
        M <- min(ris.valori)
        p.naiv <- pnorm(M, lower = TRUE)
    }
    else if (alternative == "greater") {
        M <- max(ris.valori)
        p.naiv <- pnorm(M, lower = FALSE)
    }
    else {
        M <- max(abs(ris.valori))
        p.naiv <- pnorm(M, lower = FALSE)
        onesided <- FALSE
    }
    approxx <- V * exp(-(M^2)/2)/sqrt(8 * pi)
    p.adj <- p.naiv + approxx
    p.adj <- ifelse(onesided, 1, 2) * p.adj
    out <- list(method = "Davies' test for a change in the slope",
        data.name=paste("Model = ",ogg$family$family,", link =", ogg$family$link,
        "\nformula =", as.expression(formula(ogg)),
#        data.name = paste(as.expression(ogg$call),
        "\nsegmented variable =", term),
        statistic = c("`Best' at" = valori[which.max(ris.valori)]),
        parameter = c(n.points = length(valori)), p.value = min(p.adj,1),
        alternative = alternative)
    class(out) <- "htest"
    return(out)
}


