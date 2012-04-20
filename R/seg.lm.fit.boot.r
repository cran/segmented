seg.lm.fit.boot<-function(y, XREG, Z, PSI, w, offs, opz, n.boot=10, size.boot=NULL, jt=FALSE,
    nonParam=TRUE, random=FALSE){
#random se TRUE prende valori random quando è errore: comunque devi modificare qualcosa (magari con it.max)
#     per fare restituire la dev in corrispondenza del punto psi-random
#nonParm. se TRUE implemneta il case resampling. Quello semiparam dipende dal non-errore di
#----------------------------------
#  sum.of.squares<-function(obj.seg){
#      #computes the "correct" SumOfSquares from a segmented" fit
#      b<-obj.seg$obj$coef
#      X<-qr.X(obj.seg$obj$qr) #X<-model.matrix(obj.seg)
#      X<-X[,!is.na(b)]
#      b<-b[!is.na(b)]
#      rev.b<-rev(b)
#      rev.b[1:length(obj.seg$psi)]<-0
#      b<-rev(rev.b)
#      new.fitted<-drop(X%*%b)
#      new.res<- obj.seg$obj$residuals + obj.seg$obj$fitted - new.fitted
#      ss<-sum(new.res^2)
#      ss
#      }
#--------
      n<-length(y)
      o0<-try(seg.lm.fit(y, XREG, Z, PSI, w, offs, opz), silent=TRUE)
      if(is.list(o0)){
        est.psi00<-est.psi0<-o0$psi
        ss00<-o0$SumSquares.no.gap
        } else {
        est.psi00<-est.psi0<-apply(PSI,2,mean)
        ss00<-opz$dev0
        }

      all.est.psi.boot<-all.selected.psi<-all.est.psi<-matrix(, nrow=n.boot, ncol=length(est.psi0))
      all.ss<-all.selected.ss<-rep(NA, n.boot)
      if(is.null(size.boot)) size.boot<-n

#      na<- ,,apply(...,2,function(x)mean(is.na(x)))

      Z.orig<-Z
      for(k in seq(n.boot)){
          id<-sample(n, size=size.boot, replace=TRUE)

          PSI <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
          if(jt) Z<-apply(Z.orig,2,jitter)
          if(nonParam){
              id<-sample(n, size=size.boot, replace=TRUE)
              o.boot<-try(seg.lm.fit(y[id], XREG[id,,drop=FALSE], Z[id,,drop=FALSE], PSI[id,,drop=FALSE],
                w[id], offs[id], opz), silent=TRUE)
          } else {
              yy<-o0$fitted+sample(residuals(o0),size=n, replace=TRUE)
              o.boot<-try(seg.lm.fit(yy, XREG, Z.orig, PSI, weights, offs, opz), silent=TRUE)
          }
          if(is.list(o.boot)){
            all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi

            PSI <- matrix(rep(est.psi.boot, rep(nrow(Z), length(est.psi.boot))), ncol = length(est.psi.boot))
            o<-try(seg.lm.fit(y, XREG, Z.orig, PSI, w, offs, opz), silent=TRUE)
            if(is.list(o)){
              all.est.psi[k,]<-o$psi
              all.ss[k]<-o$SumSquares.no.gap
              if(o$SumSquares.no.gap<=ifelse(is.list(o0), o0$SumSquares.no.gap, 10^12)) o0<-o
              est.psi0<-o0$psi
              all.selected.psi[k,] <- est.psi0
              all.selected.ss[k]<-o0$SumSquares.no.gap #min(c(o$SumSquares.no.gap, o0$SumSquares.no.gap))
              }
            }
          }
      all.selected.psi<-rbind(est.psi00,all.selected.psi)
      all.selected.ss<-c(ss00, all.selected.ss)

      #o<-list(all.psi.boot=all.est.psi.boot, all.psi=all.est.psi, all.selected.psi,all.ss,all.selected.ss)
      ris<-list(all.psi=drop(all.selected.psi),all.ss=all.selected.ss)
      o0$boot.restart<-ris
      return(o0)
      }

