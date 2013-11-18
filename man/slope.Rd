\name{slope}
\alias{slope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Slope estimates from segmented relationships }
\description{
  Computes the slopes of each `segmented' relationship in the fitted model.
}
\usage{
slope(ogg, parm, conf.level = 0.95, rev.sgn=FALSE, 
    var.diff=FALSE, APC=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ogg}{ an object of class "segmented", returned by any \code{segmented} method. }
  \item{parm}{ the segmented variable whose slopes have to be computed. If missing all the segmented variables are considered. }
  \item{conf.level}{ the confidence level required. }
  \item{rev.sgn}{vector of logicals. The length should be equal to the length of \code{parm}, but it is recycled 
  otherwise. When \code{TRUE} it is assumed that the current \code{parm} is `minus' the actual segmented variable,
    therefore the sign is reversed before printing. This is useful when a null-constraint has been set on the last slope.}
  \item{var.diff}{logical. If \code{var.diff=TRUE} and there is a single segmented variable, the computed standard errors 
    are based on a sandwich-type formula of the covariance matrix. See Details in \code{\link{summary.segmented}}.} 
  \item{APC}{logical. If \code{APC=TRUE} the `annual percent changes', i.e. \eqn{100\times(\exp(\beta)-1)}{100*(exp(b)-1)}, 
  are computed for each interval (\eqn{\beta}{b} is the slope). Only point estimates and confidence intervals are returned.
  } 
}
\details{
  To fit broken-line relationships, \code{segmented} uses a parameterization whose coefficients are not 
  the slopes. Therefore given an object \code{"segmented"}, \code{slope} computes point estimates,
  standard errors, t-values and confidence intervals of the slopes of each segmented relationship in the fitted model.
}
\value{
  \code{slope} returns a list of matrices. Each matrix represents a segmented relationship and its number of rows equal 
  to the number of segments, while five columns summarize the results.
}
\references{Muggeo, V.M.R. (2003) Estimating regression models with unknown break-points. 
    \emph{Statistics in Medicine} \bold{22}, 3055--3071.
    }
\author{Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} 
}
\note{The returned summary is based on limiting Gaussian distribution for the model parameters involved 
    in the computations. Sometimes, even with large sample sizes such approximations are questionable 
    (e.g., with small difference-in-slope parameters) and the results returned by \code{slope} 
    might be unreliable. Therefore is responsability of the user to gauge the applicability of such asymptotic 
    approximations. Anyway, the t values may be not assumed for testing purposes 
    and they should be used just as guidelines to assess the estimate uncertainty.
    }
\seealso{See also \code{\link{davies.test}} to test for a nonzero differece-in-slope parameter.
}
\examples{
set.seed(16)
x<-1:100
y<-2+1.5*pmax(x-35,0)-1.5*pmax(x-70,0)+rnorm(100,0,3)
out<-glm(y~1)
out.seg<-segmented(out,seg.Z=~x,psi=list(x=c(20,80)))
## the slopes of the three segments....
slope(out.seg)
rm(x,y,out,out.seg)
#
## an heteroscedastic example..
set.seed(123)
n<-100
x<-1:n/n
y<- -x+1.5*pmax(x-.5,0)+rnorm(n,0,1)*ifelse(x<=.5,.4,.1)
o<-lm(y~x)
oseg<-segmented(o,seg.Z=~x,psi=.6)
slope(oseg)
slope(oseg,var.diff=TRUE) #better CI
}
\keyword{ regression }
\keyword{ htest }
