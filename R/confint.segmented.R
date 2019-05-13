`confint.segmented` <- function(object, parm, level=0.95, method=c("delta", "score", "gradient"), rev.sgn=FALSE, var.diff=FALSE, digits=max(4, getOption("digits") - 1), ...){
#...: argomenti da passare solo a confintSegIS. Questi sono "h", "d.h", "bw" (bw="(1/n)^(1/2)"), nvalues, msgWarn o useSeg.
        method<-match.arg(method)
        cls<-class(object)
        if(length(cls)==1) cls<-c(cls, cls)
        if(method%in%c("score", "gradient") && !all(cls[1:2]==c("segmented","lm"))) stop("Score- or Gradient-based CI only work with segmented lm models") 
#=======================================================================================================
#========== metodo Delta
#=======================================================================================================
confintSegDelta<- function(object, parm, level=0.95, rev.sgn=FALSE, var.diff=FALSE, ...){
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
#--        
#        if(!"segmented"%in%class(object)) stop("A segmented model is needed")
        if(var.diff && length(object$nameUV$Z)>1) {
            var.diff<-FALSE
            warning(" 'var.diff' set to FALSE with multiple segmented variables", call.=FALSE)
            }
        #nomi delle variabili segmented:
        if(missing(parm)) {
          nomeZ<- object$nameUV$Z
          if(length(rev.sgn)==1) rev.sgn<-rep(rev.sgn,length(nomeZ))
          } else {
                if(! all(parm %in% object$nameUV$Z)) {stop("invalid 'parm' name", call.=FALSE)}
                  else {nomeZ<-parm}
                }
        if(length(nomeZ)>1) {
                warning("There are multiple segmented terms. The first is taken", call.=FALSE, immediate. = TRUE)  
                nomeZ<-nomeZ[1]
        }
        
        
        if(length(rev.sgn)!=length(nomeZ)) rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        rr<-list()
        z<-if("lm"%in%class(object)) abs(qt((1-level)/2,df=object$df.residual)) else abs(qnorm((1-level)/2))
        for(i in 1:length(nomeZ)){ #per ogni variabile segmented `parm' (tutte o selezionata)..
            #nomi.U<-grep(paste("\\.",nomeZ[i],"$",sep=""),object$nameUV$U,value=TRUE)
            #nomi.V<-grep(paste("\\.",nomeZ[i],"$",sep=""),object$nameUV$V,value=TRUE)
            nomi.U<- object$nameUV$U[f.U(object$nameUV$U, nomeZ[i])]
            nomi.V<- object$nameUV$V[f.U(object$nameUV$V, nomeZ[i])]
            m<-matrix(,length(nomi.U),3)
            colnames(m)<-c("Est.",paste("CI","(",level*100,"%",")",c(".low",".up"),sep=""))
            for(j in 1:length(nomi.U)){ #per ogni psi della stessa variabile segmented..
                    sel<-c(nomi.V[j],nomi.U[j])
                    V<-vcov(object,var.diff=var.diff)[sel,sel] #questa e' vcov di (psi,U)
                    b<-coef(object)[sel[2]] #diff-Slope
                    th<-c(b,1)
                    orig.coef<-drop(diag(th)%*%coef(object)[sel]) #sono i (gamma,beta) th*coef(ogg)[sel]
                    gammma<-orig.coef[1]
                    est.psi<-object$psi[sel[1],2]
                    V<-diag(th)%*%V%*%diag(th) #2x2 vcov() di gamma e beta
                    se.psi<-sqrt((V[1,1]+V[2,2]*(gammma/b)^2-2*V[1,2]*(gammma/b))/b^2)
                    r<-c(est.psi, est.psi-z*se.psi, est.psi+z*se.psi)
                    if(rev.sgn[i]) r<-c(-r[1],rev(-r[2:3]))
                    m[j,]<-r
                    } #end loop j (ogni psi della stessa variabile segmented)
            #CONTROLLA QUESTO:..sarebbe piu' bello
            m<-m[order(m[,1]),,drop=FALSE]
            rownames(m)<-nomi.V
            #if(nrow(m)==1) rownames(m)<-"" else m<-m[order(m[,1]),]
            if(rev.sgn[i]) {
                #m<-m[nrow(m):1,]
                rownames(m)<-rev(rownames(m))
                }
            rr[[length(rr)+1]]<- m #signif(m,digits)
            } #end loop i (ogni variabile segmented)
        names(rr)<-nomeZ
        return(rr[[1]])
          } #end_function
#=======================================================================================================
#========== metodo Score
#=======================================================================================================
confintSegIS<-function(obj, parm, d.h=1.5, h=2.5, conf.level=level, ...){
        #wrapper per ci.IS().. 
        #d.h: incremento di h..
        #se h o d.h sono negativi, tutto il range
        #==========================================================================
        #==========================================================================
        #==========================================================================
        ci.IS <- function(obj.seg, nomeZ, nomeUj, stat = c("score", "gradient"), transf=FALSE, h = -1, sigma, conf.level = 0.95, use.z = FALSE,  
                          is = TRUE, fit.is = TRUE, var.is=TRUE, bw=NULL, smooth = 0, msgWarn = FALSE, n.values = 50, 
                          altro = FALSE, cadj = FALSE, plot = FALSE, add=FALSE, agg=FALSE, raw=FALSE, useSeg=FALSE) {
                #smooth: se 0, i valori decrescenti dello IS score vengono eliminati; porta ad una curva U troppo ripida e quindi IC troppo stretti..
                #        se 2, B-spline con vincoli di monot e di "passaggio da est.psi"
                #useSeg, se TRUE (e se smooth>0) viene applicato segmented per selezionare solo i rami con pendenza negativa
                #   dovrebbe essere usato con smooth>0 e se h=-1 (all.range=TRUE)
                #transf: funziona solo con grad
                #obj.seg: oggetto restituito da segmented 
                #h: costante per definire il range dei valori di riferimento. Should be >1.
                #   Se NULL viene considerato l'intervallo 'est.psi +/- se*(zalpha*1.5) dove zalpha ? il quantile che dipende da conf.level
                #   Se qualche negativo, viene considerato il range della x dal quantile 0.02 a quello 0.98.
                #   Se >0  il range e' est.psi +/- h* zalpha * se.psi 
                # sigma se mancante viene assunta la stima presa dall'oggetto obj.seg.. 
                # use.z: se TRUE i quantili della z, otherwise la t_{n-p} 
                # stat: which statistic use 
                # agg if TRUE, and plot=TRUE and est.psi!= dalla radice che annulla lo IS score, allora l'IC ? shiftato..
                # is, fit.is, var.is: logical, induced smoothing? 
                # plot: la linea nera e' lo score originale (if raw=TRUE)
                #       la linea rossa e' lo score IS 
                #       le linea verde e' lo IS score con i pezzi decrescenti eliminati
                #       se useSeg=T aggiunge una linea segmented.. 
                #
                #          
                # conf.level: confidence levels can be vector
                # fit.is: i fitted del modello nullo provengono da un modello in cui (x-psi)_+ ?
                #         sostituito dall'approx smooth?
                # bw: the bandwidth in the kernel.. If NULL the SE(\hat\psi) is used, otherwise use a string, something like "1/n" or "sqrt(1/n)"
                # cadj: se TRUE l'approx di Ca.... che fa riferimentimento ad una Normale 
                #
                #==========================================================================
                #==========================================================================
                #==========================================================================
                u.psiX <- function(psi, sigma, x, y, XREG = NULL, scale = FALSE, est.psi = NULL, interc = FALSE, 
                                   pow = c(1, 1), lag = 0, robust = FALSE, GS = FALSE, is = FALSE, se.psi, var.is = TRUE, which.return = 3,
                                   fit.is = FALSE, altro = FALSE, cadj = FALSE, transf=FALSE) {
                        # Restituisce score e/o var, e/o score stand. (vedi 'which.return') Inoltre se robust=TRUE calcola la
                        # var robusta est.psi: o NULL oppure uno scalare con attributi 'b' e 'fitted' se lag>0 allora la
                        # variabile V viene modificata nell'intorno di psi. Valori di pow diversi da uno sono ignorati quando
                        # lag>0 pow: due potenze dei termini (x-psi)_+ e I(x>psi) se GS=TRUE calcola la statistica GS.
                        # richiede 'est.psi', e 'scale' ? ignorato which.return. 3 means the scaled score, 1= the unscaled
                        # score, 2=the sqrt(variance) (see the last row) 
                        # is: se TRUE lo smoothing indotto al num 
                        # var.is: se TRUE lo smooth indotto viene usato anche per il denom (ovvero per la var dello score)
                        # U.is: se TRUE (provided that is=TRUE) the design matrix includes (x-psi)*pnorm((x-psi)/se) rather than pmax(x-psi,0)
                        #altro: se TRUE (and fit.is=TRUE), U.psi = (x-psi)*pnorm((x-psi)/se) + h*dnorm((x-psi)/h)
                        #--------------------------------------------
                        
                        varUpsi.fn <- function(X, sigma = 1, r = NULL) {
                                #X: the design matrix. The 1st column corresponds to psi 
                                #r: the residual vector. If NULL the usual model-based (rather than robust) variance is returned. INF<- if(length(sigma)==1)
                                # (sigma^2)*crossprod(X) else crossprod(X,diag(sigma^2))%*%X
                                INF <- crossprod(X)/(sigma^2)
                                if (is.null(r)) {
                                        vv <- INF[1, 1] - (INF[1, -1] %*% solve(INF[-1, -1], INF[-1, 1]))
                                } else {
                                        u <- X * r/(sigma^2)
                                        V <- crossprod(u)  #nrow(X)*var(u)
                                        I22 <- solve(INF[-1, -1])
                                        vv <- V[1, 1] - INF[1, -1] %*% I22 %*% V[1, -1] - V[1, -1] %*% I22 %*% INF[-1, 1] + INF[1,
                                                                                                                                -1] %*% I22 %*% V[-1, -1] %*% I22 %*% INF[-1, 1]
                                }
                                return(vv)
                        }
                        # f.f<-function(x,psi,l=0){ x1<-1*I(x>psi) id<-which(x1>=1)[1] id.change <-
                        # max(1,(id-l)):min(length(x),(id+l)) val<-((1/(2*l+1))*( 1:(2*l+1)))[1:length(id.change)]
                        # #if(length(id.change)!=length(val)) return x1[id.change]<-val x1<- -x1 x1 }
                        dpmax <- function(x, y, pow = 1) {
                                # derivata prima di pmax; se pow=1 ? -I(x>psi)
                                if (pow == 1)
                                        ifelse(x > y, -1, 0) else -pow * pmax(x - y, 0)^(pow - 1)
                        }
                        if (cadj && which.return != 3)
                                stop("cadj=TRUE can return only the studentized score")
                        if (is && missing(se.psi))
                                stop("is=TRUE needs se.psi")
                        if (interc) XREG <- cbind(rep(1, length(y)), XREG)
                        
                        if(fit.is) {
                                XX<- if(altro) cbind((x-psi)*pnorm((x - psi)/se.psi)+se.psi*dnorm((x-psi)/se.psi), XREG) else cbind((x-psi)*pnorm((x-psi)/se.psi), XREG)
                                o <- lm.fit(x = XX, y = y)
                                #o <- lm.fit(x = cbind(XREG, (x - psi) * pnorm((x - psi)/se.psi)), y = y)
                        } else {
                                XX<- cbind(pmax(x - psi, 0)^pow[1], XREG)
                                o <- lm.fit(x = XX, y = y)  
                                #o <- lm.fit(x = cbind(XREG, pmax(x - psi, 0)), y = y)  #o<-lm(y~0+XREG+pmax(x-psi,0))
                        }
                        
                        #b <- o$coef[length(o$coef)]
                        b <- o$coef[1]
                        mu <- o$fitted.values
                        n <- length(mu)
                        #  if (cadj) sigma <- sqrt(sum(o$residuals^2)/(n - sum(!is.na(o$coef)) - 1))
                        #  V <- if (lag == 0) dpmax(x, psi, pow = pow[2]) else f.f(x, psi, lag)  #V <- rowMeans(sapply(x, function(xx){-I(x>xx)}))
                        V<-NULL #serve per il check..
                        if (GS) {
                                if (is.null(est.psi)) stop("'GS=TRUE' needs 'est.psi'")
                                gs <- b * (sum((y - mu) * V)/(sigma^2)) * (est.psi - psi)
                                gs <- sqrt(pmax(gs, 0)) * sign(est.psi - psi)
                                return(gs)
                        }
                        if(is){
                                r<- -b*sum(((y-mu)*pnorm((x - psi)/se.psi)))/sigma^2       
                                XX<- if(var.is) cbind(-b*pnorm((x - psi)/se.psi), XX) else cbind(-b*I(x > psi), XX)
                        } else {
                                r<- -b*sum((y-mu)*I(x > psi))/sigma^2
                                XX<- cbind(-b*I(x > psi), XX)
                        }
                        
                        #XX <- if (is) cbind(-b * pnorm((x - psi)/se.psi), (x - psi)*pnorm((x - psi)/se.psi), XREG) else cbind(b * V, pmax(x - psi, 0)^pow[1], XREG)
                        #r <- drop(crossprod(XX, y - mu))/sigma^2
                        #if (is && altro) r[1] <- r[1] + (b^2) * se.psi * sum(dnorm((x - psi)/se.psi))/sigma^2
                        #if (!var.is) XX <- cbind(b * V, pmax(x - psi, 0)^pow[1], XREG)
                        if (scale) {
                                if (!is.null(est.psi)) {
                                        # questo e' se devi usare l'inf osservata. Cmq visto che dipende da est.psi e non psi, se scale=TRUE
                                        # sarebbe inutile calcolarla ogni volta..
                                        mu <- attr(est.psi, "fitted")
                                        est.b <- attr(est.psi, "b")
                                        est.psi <- as.numeric(est.psi)
                                        #V <- if (lag == 0) dpmax(x, est.psi, pow = pow[2]) else f.f(x, est.psi, lag)  #V <- rowMeans(sapply(x, function(xx){-I(x>xx)}))
                                        #XX <- cbind(est.b * V, pmax(x - psi, 0)^pow[1], XREG)
                                        if(is){
                                                XX<- if(var.is) cbind(-est.b*pnorm((x - est.psi)/se.psi), XX[,-1]) else cbind(-est.b*I(x > est.psi), XX[,-1])
                                        } else {
                                                XX<- cbind(-est.b*I(x > est.psi), XX[,-1])
                                        }
                                }
                                # INF<- if(length(sigma)==1) (sigma^2)*crossprod(XX) else crossprod(XX,diag(sigma^2))%*%XX
                                # v.Upsi<-INF[1,1]-(INF[1,-1] %*% solve(INF[-1,-1],INF[-1,1]))
                                rr <- if (robust) (y - mu) else NULL
                                v.Upsi <- try(varUpsi.fn(XX, sigma, r = rr), silent = TRUE)
                                if (!is.numeric(v.Upsi))
                                        return(NA)
                                if (v.Upsi <= 0)
                                        return(NA)
                                # r<-r[1]/sqrt(v.Upsi)
                        }
                        names(r) <- NULL
                        #r <- c(r[1], v.Upsi, r[1]/sqrt(max(v.Upsi, 0)))
                        r <- c(r, v.Upsi, r/sqrt(max(v.Upsi, 0)))
                        r <- r[which.return]
                        if (cadj)
                                r <- sign(r) * sqrt((r^2) * (1 - (3 - (r^2))/(2 * n)))
                        r
                }
                
                # per disegnare devi vettorizzare
                u.psiXV <- Vectorize(u.psiX, vectorize.args = "psi", USE.NAMES = FALSE)
                
                #==========================================================================
                gs.fn <- function(x, y, estpsi, sigma2, psivalue, pow = c(1,1), adj = 1, is = FALSE,
                                  sepsi, XREG = NULL, fit.is = FALSE, altro = FALSE, transf=FALSE) {
                        # calcola la statist gradiente 
                        #x,y i dati; estpsi la stima di psi 
                        #a: la costante per lisciare I(x>psi)-> aI(x>psi)^{a-1} (ignorata se is=TRUE)
                        #  
                        # is: se TRUE calcola la GS usando lo score 'naturally smoothed' 
                        #adj. Se 0 non fa alcuna modifica e cosi' potrebbe risultare non-positiva.  Se 1 e 2 vedi i codici all'interno
                        logitDeriv<-function(kappa) exp(kappa)*diff(intv)/((1+exp(kappa))^2)
                        logit<-function(psi) log((psi-min(intv))/(max(intv)-psi))
                        logitInv<-function(kappa) (min(intv)+max(intv)*exp(kappa))/(1+exp(kappa))
                        intv<-quantile(x, probs=c(.02,.98),names=FALSE)
                        if (is && missing(sepsi))
                                stop("SE(psi) is requested when is=TRUE")
                        k <- length(psivalue)
                        r <- vector(length = k)
                        for (i in 1:k) {
                                psii <- psivalue[i]
                                #prima dell'aggiunta di altro..'
                                #    if (fit.is) {
                                #      X <- cbind(1, x, (x - psii) * pnorm((x - psii)/sepsi), XREG)
                                #    } else {
                                #      X <- cbind(1, x, pmax(x - psii, 0), XREG)
                                #    }
                                
                                if(fit.is) {
                                        X<- if(altro) cbind(1,x, (x-psii)*pnorm((x - psii)/sepsi)+sepsi*dnorm((x-psii)/sepsi), XREG) else cbind(1,x,(x-psii)*pnorm((x-psii)/sepsi), XREG)
                                } else {
                                        X<- cbind(1,x,pmax(x - psii, 0)^pow[1], XREG)
                                }
                                
                                o <- lm.fit(y = y, x = X)
                                b <- o$coef[3]
                                if (is) {
                                        v <- pnorm((x - psii)/sepsi)
                                } else {
                                        v <- if (pow[2] == 1) I(x > psii) else pow[2] * pmax(x - psii, 0)^(pow[2] - 1)
                                }
                                if(transf) v<-v * logitDeriv(logit(psii))
                                r[i] <- -(b/sigma2) * sum((y - o$fitted) * v) 
                                r[i] <- if(!transf) r[i]*(estpsi - psii) else r[i]*(logit(estpsi) - logit(psii))
                                if (altro && fit.is)
                                        r[i] <- r[i] + (estpsi - psii) * ((b * sepsi * sum(dnorm((x - psii)/sepsi))) * (b/sigma2))
                        }
                        if (adj > 0) {
                                r<- if (adj == 1) pmax(r, 0) else abs(r)
                        }
                        
                        if(transf) psivalue<-logit(psivalue)
                        segni<-if(transf) sign(logit(estpsi) - psivalue) else sign(estpsi - psivalue) 
                        #plot(psivalue, r, type="o")
                        r <- cbind(psi = psivalue, gs.Chi = r, gs.Norm = sqrt(r) * segni )
                        r
                }
                #==========================================================================
                monotSmooth <- function(xx, yy, hat.psi, k = 20, w = 0) {
                        # xx: esplicativa yy: yy la risposta hat.psi: la stima del psi k: se ? uno scalare allora il rango
                        # della base, altrimenti i nodi.. w: l'esponente per costruire il vettore dei pesi (per dare pi? peso
                        # 'localmente')
                        #-------------------
                        bspline <- function(x, ndx, xlr = NULL, knots, deg = 3, deriv = 0) {
                                # x: vettore di dati xlr: il vettore di c(xl,xr) ndx: n.intervalli in cui dividere il range deg: il
                                # grado della spline
                                #require(splines)
                                if (missing(knots)) {
                                        if (is.null(xlr)) {
                                                xl <- min(x) - 0.01 * diff(range(x))
                                                xr <- max(x) + 0.01 * diff(range(x))
                                        } else {
                                                if (length(xlr) != 2)
                                                        stop("quando fornito, xlr deve avere due componenti")
                                                xl <- xlr[1]
                                                xr <- xlr[2]
                                        }
                                        dx <- (xr - xl)/ndx
                                        knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
                                }
                                B <- splineDesign(knots, x, ord = deg + 1, derivs = rep(deriv, length(x)))
                                # B<-spline.des(knots,x,bdeg+1,0*x) #$design
                                r <- list(B = B, degree = deg, knots = knots)  #, dx=dx, nterm=ndx)
                                r  #the B-spline base matrix
                        }  #end_fn
                        #---------
                        if (length(k) == 1)
                                r <- bspline(xx, ndx = k) else r <- bspline(xx, knots = k)
                                B <- r$B
                                knots <- r$knots
                                degree <- r$degree
                                
                                D1 <- diff(diag(ncol(B)), diff = 1)
                                d <- drop(solve(crossprod(B), crossprod(B, yy)))
                                # calcola monotone splines. La pen si riferisce solo alle diff dei coef della base!!
                                
                                # rx <- range(xx) nterm <- round(nterm) dx <- (rx[2] - rx[1])/nterm knots <- c(rx[1] + dx *
                                # ((-degree):(nterm - 1)), rx[2] + dx * (0:degree))
                                
                                B0 <- spline.des(knots, c(min(xx), hat.psi, max(xx)), degree + 1)$design
                                P <- tcrossprod(B0[2, ]) * 10^12
                                
                                e <- rep(1, length(d))
                                ww <- (1/(abs(xx - hat.psi) + diff(range(xx))/100))^w
                                it <- 0
                                while (!isTRUE(all.equal(e, rep(0, length(e))))) {
                                        v <- 1 * I(diff(d) > 0)
                                        E <- (10^12) * crossprod(D1 * sqrt(v))  #t(D1) %*%diag(v)%*%D1 #
                                        d.old <- d  #a.new
                                        M <- crossprod(B * sqrt(ww)) + E + P  #t(B)%*% B + E + P
                                        d <- drop(solve(M+.001*diag(ncol(M)), crossprod(B, ww * yy)))  #d <- drop(solve(M, t(B)%*% yy))
                                        e <- d - d.old
                                        it <- it + 1
                                        if (it >= 20)
                                                break
                                }  #end_while
                                fit <- drop(B %*% d)
                                return(fit)
                }
                #==========================================================================
                miop<-function(x,y,xs=x, ys=y, h=FALSE,v=FALSE, only.lines=FALSE,
                               top=TRUE, right=TRUE, col.h=grey(.6), col.v=col.h,...){
                        #disegna il calssico plot(x,y,..) e poi aggiunge le proiezioni orizzontali e/o verticali
                        #x, y : vettori per cui disegnare il grafico
                        #xs, ys: punti rispetto a cui disegnare le proiezioni (default a tutti)
                        #h, v: disegnare le linee horizontal and vertical?
                        #top: le linee v riportarle verso l'alto (TRUE) o il basso?
                        #right: le linee horiz riportarle verso destra (TRUE) o sinistra?
                        #only.lines: se TRUE disegna (aggiungendo in un plot *esistente*) solo le "proiezioni" (linee "v" e "h")
                        if(only.lines) h<-v<-TRUE
                        if(!only.lines) plot(x,y,type="l",...)
                        #  col.h<-col.v<-1:length(xs)
                        if(v){
                                y0<- if(top) par()$usr[4] else par()$usr[3]
                                segments(xs, y0, xs,ys, col=col.v, lty=3)
                        }
                        if(h){
                                x0<-if(right) par()$usr[2] else  par()$usr[1]
                                segments(xs,ys, x0,ys,col=col.h, lty=3, lwd=1.2)
                        }
                        invisible(NULL)
                }
                #==========================================================================
                f.Left<-function(x,y){
                        yy<-rev(y)
                        xx<-rev(x)
                        idList<-NULL
                        while(any(diff(yy)<0)){
                                id<-which(diff(yy)<0)[1]
                                idList[length(idList)+1]<- id+1
                                yy<-yy[-(id+1)]
                                xx<-xx[-(id+1)]
                        }
                        r<-cbind(xx,yy)
                        r
                }
                #==========================================================================
                f.Right<-function(x,y){
                        #elimina i valori che violano la monotonic
                        xx<-x
                        yy<-y
                        idList<-NULL
                        while(any(diff(yy)>0)){
                                id<-which(diff(yy)>0)[1]
                                idList[length(idList)+1]<- id+1
                                yy<-yy[-(id+1)]
                                xx<-xx[-(id+1)]
                        }
                        r<-cbind(xx,yy)
                        r
                }
                #==========================================================================
                #==========================================================================
                #==========================================================================
                
                stat <- match.arg(stat)
                if (missing(sigma)) sigma <- summary.lm(obj.seg)$sigma
                if (cadj) use.z = TRUE
                zalpha <- if (use.z) -qnorm((1 - conf.level)/2) else -qt((1 - conf.level)/2, df = obj.seg$df.residual)
                
                if(!is.numeric(h)) stop(" 'h' should be numeric")
                if(sign(h)>=0) h<-abs(h[1])
                
                Y <- obj.seg$model[, 1]  #la risposta
                X <- obj.seg$model[, nomeZ]
                formula.lin<- update.formula(formula(obj.seg), paste(".~.", paste("-",paste(obj.seg$nameUV$V,collapse =  "-")))) #remove *all* V variables
                formula.lin<- update.formula(formula.lin, paste(".~.-", nomeUj))
                #formula.lin <- update.formula(formula(obj.seg), paste(".~.", paste("-",paste(c(obj.seg$nameUV$U,obj.seg$nameUV$V),collapse =  "-"))))
                XREG <- model.matrix(formula.lin, data = obj.seg$model)
                if (ncol(XREG) == 0) XREG <- NULL
                                             
                nomePsij<-sub("U","psi", nomeUj)
                est.psi <- obj.seg$psi[nomePsij, "Est."]
                se.psi <- obj.seg$psi[nomePsij, "St.Err"]
                                             
                if (any(h < 0)) {
                     all.range <- TRUE
                     valori <- seq(quantile(X,probs=.05, names=FALSE), quantile(X,probs=.95, names=FALSE), l = n.values)
                     } else {
                     all.range <- FALSE
                     valori <- seq(max(quantile(X,probs=.05, names=FALSE), est.psi - h * se.psi), 
                                min(quantile(X,probs=.95, names=FALSE), est.psi + h * se.psi), l = n.values)
                     }
                n <- length(Y)
                min.X <- min(X)
                max.X <- max(X)
                if(!is.null(bw)) se.psi<-eval(parse(text=bw))
                if (stat == "score") {
                        U.valori <- u.psiXV(psi = valori, sigma = sigma, x = X, y = Y, XREG = XREG, is = is, se.psi = se.psi,
                                scale = TRUE, pow = c(1, 1), fit.is = fit.is, altro = altro,  cadj = cadj, var.is=var.is, transf=transf)
                        statlab<-"Score statistic"
                        if(plot && raw)  U.raw <- u.psiXV(valori, sigma, X, Y, XREG, is=FALSE, scale=TRUE, pow = c(1, 1), fit.is=FALSE, altro =altro, cadj = cadj, var.is=FALSE, transf=transf)  
                        } else {
                        U.valori <- gs.fn(X, Y, est.psi, sigma^2, valori, is = is, sepsi = se.psi, XREG = XREG, 
                                fit.is = fit.is, altro = altro, transf=transf, pow=c(1,1))[, 3]
                        statlab<-"Gradient statistic"
                        if(plot && raw)  U.raw <- gs.fn(X, Y, est.psi, sigma^2, valori, is=FALSE, XREG=XREG, fit.is=FALSE, altro=altro, transf=transf)[,3]  
                                }
                        
                if(any(is.na(U.valori))) { #stop("NA in the statistic values")	
                        warning("removing NA in the statistic values")
                        valori<-valori[!is.na(U.valori)]	
                        U.valori<-U.valori[!is.na(U.valori)]
                        }
                                             
                logit<-function(psi) log((psi-min(intv))/(max(intv)-psi))
                logitInv<-function(kappa) (min(intv)+max(intv)*exp(kappa))/(1+exp(kappa))
                intv<-quantile(X, probs=c(.02,.98),names=FALSE)

                if (stat == "gradient" && transf) {
                        est.psi<- logit(est.psi)
                        valori<- logit(valori)
                        x.lab<- "kappa"
                        }                
                                             
                if(plot && !add) {
                        x.lab<-"psi"
                        if(raw) {
                                plot(valori, U.raw, xlab=x.lab, ylab=statlab, type="l")
                                points(valori, U.valori, xlab=x.lab, ylab=statlab, type="l", col=2) 
                                } else {
                                plot(valori, U.valori, xlab=x.lab, ylab=statlab, type="l", col=2)
                                }
                                abline(h=0, lty=3)
                                segments(est.psi,0, est.psi, -20, lty=2)
                                }
                                             
                if(prod(range(U.valori))>=0) stop("the signs of stat at extremes are not discordant, increase 'h' o set 'h=-1' ")    
                
                if(smooth==0){
                        #rimuovi i pezzi di U.valori decrescenti..
                        ####left
                        valoriLeft<-valori[valori<=est.psi]  #valori[U.valori>=0]
                        UvaloriLeft<-U.valori[valori<=est.psi] #U.valori[U.valori>=0]
                        vLeft<-f.Left(valoriLeft,UvaloriLeft) #rendi monotona la curva..
                        valori.ok<-vLeft[,1]
                        Uvalori.ok<-vLeft[,2]
                        f.interpL <- splinefun(Uvalori.ok, valori.ok, method="mono")
                        ####right
                        valoriRight<-valori[valori>=est.psi]  #valori[U.valori<0]
                        UvaloriRight<-U.valori[valori>=est.psi] #U.valori[U.valori<0]
                        vRight<-f.Right(valoriRight,UvaloriRight)
                        valori.ok<-vRight[,1]
                        Uvalori.ok<-vRight[,2]
                        f.interpR <- splinefun(Uvalori.ok, valori.ok, method="mono")
                                } else { #if smooth>0
                        if(useSeg){
                           oseg<-try(suppressWarnings(segmented(lm(U.valori~valori), ~valori, psi=quantile(valori, c(.25,.75),names=FALSE), 
                                control=seg.control(n.boot=0,stop.if.error = F))),silent=TRUE)
                           #seg.lm.fit.boot(U.valori, XREG, Z, PSI, w, offs, opz)
                           if(class(oseg)[1]=="try-error"){
                                oseg<-try(suppressWarnings(segmented(lm(U.valori~valori), ~valori, psi=quantile(valori, .5,names=FALSE), 
                                        control=seg.control(n.boot=0))),silent=TRUE)
                                        }
                           if(class(oseg)[1]=="segmented"){
                                if(plot) lines(valori, oseg$fitted, lty=3, lwd=1.5)
                                soglie<-oseg$psi[,2]
                                iid<-cut(valori,c(min(valori)-1000, soglie, max(valori)+1000), labels=FALSE)
                                slopes<-cumsum(oseg$coef[2:(length(oseg$coef)-length(soglie))])
                                slopes<-rep(slopes,table(iid))
                                valori<-valori[slopes<=0]
                                U.valori<-U.valori[slopes<=0]
                                }
                           } 
                        fr<-monotSmooth(valori,U.valori,est.psi,k=7)
                        fr<- fr -(.2/diff(range(valori))) *(valori-mean(valori)) #add a small negative trend to avoid constant values in U..
                        vLeft<-cbind(valori[valori<=est.psi], fr[valori<=est.psi])
                        vRight<-cbind(valori[valori>=est.psi], fr[valori>=est.psi])
                        if(!all.range){
                                if( (min(valori)> intv[1]) && (fr[1]< max(zalpha))) return("errLeft")
                                if( (max(valori)< intv[2]) && (fr[length(fr)]> min(-zalpha))) return("errRight")
                                }
                        f.interpL<-f.interpR<-splinefun(fr,valori,"m")
                        }#end_if smooth 
                L<-f.interpL(zalpha) 
                U<-f.interpR(-zalpha)
                #browser()    
                #il valore che annulla lo IS score puo' essere differente dalla stima di segmented
                #   quindi salviamo questo "delta": gli IC potrebbero essere aggiustati con IC+delta
                delta<- est.psi-f.interpL(0)  #if(abs((f.interpL(0)-f.interpR(0))/f.interpR(0))>.001)
                                             
                if(plot){
                        if(!agg) delta<-0
                        #if(raw) plot(valori, U.raw, xlab="psi", ylab=statlab, type="l") else plot(valori, U.valori, xlab="psi", ylab=statlab, type="n")
                        lines(vLeft, col=3); lines(vRight, col=3)
                        vv<-seq(0,zalpha*1.2,l=50)
                        lines(f.interpL(vv)+delta,vv, col=grey(.8, alpha=.6), lwd=4)
                        vv<-seq(0,-zalpha*1.2,l=50)
                        lines(f.interpR(vv)+delta,vv, col=grey(.8, alpha=.6), lwd=4)
                        points(est.psi, 0, pch=19)             
                        miop(c(L,U)+delta,c(zalpha,-zalpha),only.lines=TRUE,top=FALSE, right=FALSE)
                        }
                if (stat == "gradient" && transf) {
                        L<-logitInv(L)
                        U<-logitInv(U)
                        }
                L<- pmax(L, quantile(X,probs=.02))
                U<- pmin(U,quantile(X,probs=.98))
                                             
                #r<-cbind(lower=L,upper=U)
                #rownames(r) <- paste(conf.level)
                #attr(r, "delta")<-delta
                r<-c(est.psi, L, U)
                return(r)
                } #end fn
        
        #--------------------------------------------------------------------------
        #==========================================================================
        #==========================================================================
        #==========================================================================

        if(!all(class(obj) == c("segmented","lm"))) stop("A segmented lm object is requested")
        if(missing(parm)){
                nomeZ<- parm<- obj$nameUV$Z
                } else {
                if(!all(parm %in% obj$nameUV$Z)) stop("invalid 'parm' ")
                nomeZ<-parm
                }
        if(length(parm)>1) {
                warning("There are multiple segmented terms. The first is taken", call.=FALSE, immediate. = TRUE)
                nomeZ<-parm[1]
                }
        nomiU.term<-grep(nomeZ, obj$nameUV$U, value=TRUE) #termini U per la *stessa* variabile..
        #npsi.term<- length(nomiU.term) #no. di breakpoints for the same variable.
        ra<-matrix(NA, length(nomiU.term), 3)
        rownames(ra)<- nomiU.term
        
        for(U.j in nomiU.term){
                if(any(c(d.h, h)<0)) {
                        ra[U.j,]<-ci.IS(obj, nomeZ, U.j, h=-1, conf.level=level, ...)
                }  
                d.h<-min(max(d.h, 1.5),10)
                a<-"start"
                it<-0
                while(is.character(a)){
                        a<- try(ci.IS(obj, nomeZ, U.j, h=h, conf.level=level, ...), silent=TRUE)
                        h<-h*d.h
                        it<-it+1
                        #cat(it,"\n")
                        if(it>=20) break
                        }
                ra[U.j,]<-a
                }
        colnames(ra)<-c("Est.",paste("CI","(",level*100,"%",")",c(".low",".up"),sep=""))
        rownames(ra)<-sub("U","psi", nomiU.term)
        ra
        } #end fn confintSegIS
        
#=======================================================================================================
#========== inizio funzione
#=======================================================================================================

if(method=="delta"){
        r<-confintSegDelta(object, parm, level, rev.sgn, var.diff, ...)    
        } else {
        r<-confintSegIS(object, parm, stat=method, conf.level=level, ...)       
        }
r<-signif(r,digits)
return(r)                
}
