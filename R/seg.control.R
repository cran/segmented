`seg.control` <-
function(n.boot=10, display=FALSE, tol=1e-5, it.max=30, fix.npsi=TRUE, K=10, quant=TRUE, maxit.glm=25, h=1, break.boot=5,
         size.boot=NULL, jt=FALSE, nonParam=TRUE, random=TRUE, seed=12345, fn.obj=NULL, digits=NULL, conv.psi=FALSE, alpha=.02, min.step=.0001,     
         powers=c(1,1), last=TRUE, stop.if.error=NULL, gap=FALSE, fc=.95){
      list(toll=tol,it.max=it.max,visual=display,stop.if.error=stop.if.error,
            K=K,last=last,maxit.glm=maxit.glm,h=h,n.boot=n.boot, size.boot=size.boot, gap=gap, jt=jt, break.boot=break.boot,
            nonParam=nonParam, random=random, pow=powers, seed=seed, quant=quant, fn.obj=fn.obj, digits=digits, 
            conv.psi=conv.psi, alpha=alpha, fix.npsi=fix.npsi, min.step=min.step, fc=fc)}
