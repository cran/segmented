`davies.test` <- 
function(ogg, term, k=10, alternative="two-sided"){
    if (!is.character(term))
        stop("term must be character")
    Z<-if(term%in%names(ogg$model)) {
        ogg$model[,term]} else {eval(parse(text = term))}
    qq <- quantile(Z, prob = c(0.05, 0.95), names = FALSE,na.rm=TRUE)
    sx <- qq[1]
    dx <- qq[2]
    valori <- seq(sx, dx, length = k)
    ris.valori <- NULL
    for(i in valori){
        o <- segmented(ogg, seg.Z = ~Z, psi = list(Z = i), control = seg.control(it.max = 0,
            display = FALSE))
        ris.valori[(length(ris.valori)) + 1] <- if (is.list(o)) {
            summary(o)$coef[length(coef(ogg))+1,3]
            } else {    NA    }
    } #end_for
    valori <- valori[!is.na(ris.valori)]
    ris.valori <- ris.valori[!is.na(ris.valori)]
    V <- sum(abs(diff(ris.valori)))
    onesided <- alternative == "one-sided"
    if (!onesided) ris.valori <- abs(ris.valori)
    M <- max(ris.valori)
    approxx <- V * exp(-(M^2)/2)/sqrt(8 * pi)
    p.naiv <- pnorm(-M)
    p.adj <- p.naiv + approxx
    p.adj <- ifelse(onesided, 1, 2) * p.adj
    out <- list(method = "Davies' test for a change in slope",
        data.name = paste("model =", deparse(substitute(ogg)),
            ",", "segmented variable =", term), statistic = c("maximum at" = valori[which.max(ris.valori)]),
        parameter = c(n.points = length(valori)), p.value = min(p.adj,
            1), alternative = alternative)
    class(out) <- "htest"
    return(out)
    }
