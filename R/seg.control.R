`seg.control` <-
function(n.boot=10, display=FALSE, tol=1e-5, it.max=30, fix.npsi=TRUE, K=10, quant=FALSE, maxit.glm=NULL, h=1.25, break.boot=5,
         size.boot=NULL, jt=FALSE, nonParam=TRUE, random=TRUE, seed=NULL, fn.obj=NULL, digits=NULL, alpha=NULL, 
         fc=.95, check.next=TRUE, tol.opt=NULL, fit.psi0=NULL, eta=NULL,min.nj=2){
      list(toll=tol,it.max=it.max,visual=display,stop.if.error=NULL,
            K=K,last=TRUE, maxit.glm=maxit.glm,h=h,n.boot=n.boot, size.boot=size.boot, gap=FALSE, jt=jt, break.boot=break.boot,
            nonParam=nonParam, random=random, pow=c(1,1), seed=seed, quant=quant, fn.obj=fn.obj, digits=digits, 
            conv.psi=FALSE, alpha=alpha, fix.npsi=fix.npsi, fc=fc, check.next=check.next, 
           tol.opt=tol.opt, fit.psi0=fit.psi0, eta=eta,min.n=min.nj)}

