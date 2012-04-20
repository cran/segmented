`seg.control` <-
function(toll=.0001, it.max=10, display=FALSE, stop.if.error=TRUE, K=10, last=TRUE, maxit.glm=25, h=1, 
    n.boot=10, size.boot=NULL, gap=FALSE, jt=FALSE, nonParam=TRUE, random=TRUE){
        list(toll=toll,it.max=it.max,visual=display,stop.if.error=stop.if.error,
            K=K,last=last,maxit.glm=maxit.glm,h=h,n.boot=n.boot, size.boot=size.boot, gap=gap, jt=jt,
            nonParam=nonParam, random=random)}

