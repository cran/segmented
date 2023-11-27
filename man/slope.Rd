\name{slope}
\alias{slope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Slope estimates from segmented/stepmented relationships }
\description{
  Computes the slopes of each `segmented' (or even `stepmented') relationship in the fitted model.
}
\usage{
slope(ogg, parm, conf.level = 0.95, rev.sgn=FALSE, 
    APC=FALSE, .vcov=NULL, .coef=NULL, 
    use.t=NULL, by=NULL, interc=TRUE, ..., digits = max(4, getOption("digits") - 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ogg}{ an object of class "segmented", returned by any \code{segmented} method or a list of two segmented fits to compare the estimates of corresponding slopes.}
  \item{parm}{ the segmented variable whose slopes have to be computed. If missing all the segmented variables are considered. }
  \item{conf.level}{ the confidence level required. }
  \item{rev.sgn}{vector of logicals. The length should be equal to the length of \code{parm}, but it is recycled 
  otherwise. When \code{TRUE} it is assumed that the current \code{parm} is `minus' the actual segmented variable,
    therefore the sign is reversed before printing. This is useful when a null-constraint has been set on the last slope.}
  \item{APC}{logical. If \code{APC=TRUE} the `annual percent changes', i.e. \eqn{100\times(\exp(\beta)-1)}{100*(exp(b)-1)}, 
  are computed for each interval (\eqn{\beta}{b} is the slope). Only point estimates and confidence intervals are returned.}
  \item{.vcov}{ The \emph{full} covariance matrix of estimates. If unspecified (i.e. \code{NULL}), the covariance matrix is computed internally by \code{vcov(ogg)}.}
  \item{.coef}{ The regression parameter estimates. If unspecified (i.e. \code{NULL}), it is computed internally by \code{coef(ogg)}.}
  \item{use.t}{Which quantiles should be used to compute the confidence intervals? If \code{NULL} (default) the \eqn{t}{t} distribution is used only for objects obtained by \code{segmented.lm}. }
  \item{by}{Only for \code{segmented.lme} objects. 
  It is a named list indicating covariate names and corresponding values affecting the fitted segmented relationship. 
  For instance, \code{by=list(group="2",z2=.2)}, provided that the model has been fitted by specifying \code{group} and \code{z2} in \code{x.diff} (or as interaction with the segmented variable). Note that if the provided variables or values are irrelevant for changing the slopes, a warning message is printed.   
  }
  \item{interc}{logical, only for \code{'stepmented'} fits. If \code{TRUE}, the mean levels also account for the intercept; otherwise the first level is assumed to be zero.
  }
  \item{...}{ Further arguments to be passed on to \code{vcov.segmented}, such as \code{var.diff} and \code{is}. See Details in \code{\link{vcov.segmented}} and \code{\link{summary.segmented}}.
  }
\item{digits}{controls number of digits in the returned output.}
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
\references{
    Muggeo, V.M.R. (2003) Estimating regression models with unknown break-points. 
    \emph{Statistics in Medicine} \bold{22}, 3055--3071.
    }
\author{Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} }
\note{The returned summary is based on limiting Gaussian distribution for the model parameters involved 
    in the computations. Sometimes, even with large sample sizes such approximations are questionable 
    (e.g., with small difference-in-slope parameters) and the results returned by \code{slope} 
    might be unreliable. Therefore is responsability of the user to gauge the applicability of such asymptotic 
    approximations. Anyway, the t values may be not assumed for testing purposes 
    and they should be used just as guidelines to assess the estimate uncertainty.
    }
\seealso{See also \code{\link{davies.test}} and \code{\link{pscore.test}} to test for a nonzero difference-in-slope parameter.
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
