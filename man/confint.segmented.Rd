\name{confint.segmented}
\alias{confint.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence intervals for breakpoints}
\description{
  Computes confidence intervals for the breakpoints in a fitted `segmented' model.
}
\usage{
\method{confint}{segmented}(object, parm, level=0.95, rev.sgn=FALSE, var.diff=FALSE,
        digits=max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted \code{segmented} object. }
  \item{parm}{the segmented variable of interest. If missing all the segmented variables are considered. }
  \item{level}{the confidence level required (default to 0.95).}
  \item{rev.sgn}{vector of logicals. The length should be equal to the length of \code{parm}; recycled otherwise.
   when \code{TRUE} it is assumed that the current \code{parm} is `minus' the actual segmented variable,
    therefore the sign is reversed before printing. This is useful when a null-constraint has been set on the last slope.}
  \item{var.diff}{logical. If \code{var.diff=TRUE} and there is a single segmented variable, the standard error is 
    based on sandwich-type formula of the covariance matrix. See Details in \code{\link{summary.segmented}}.}
  \item{digits}{controls the number of digits to print when printing the output. }
  \item{\dots}{additional parameters }
}
\details{
  Currently \code{confint.segmented} computes confidence limits for the breakpoints using the standard error coming from the Delta
  method for the ratio of two random variables. This value is an approximation (slightly) better than the 
  one reported in the `psi' component of the list returned by any \code{segmented} method. The resulting
  confidence intervals are based on the asymptotic Normal distribution of the breakpoint
  estimator which is reliable just for clear-cut kink relationships. See Details in \code{\link{segmented}}.
  }
\value{
    A list of matrices. Each matrix includes point estimate and confidence limits of the breakpoint(s) for each
    segmented variable in the model.
}
%\references{ }
\author{ Vito M.R. Muggeo }
%\note{ ~~further notes~~
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{segmented}} }
\examples{
set.seed(10)
x<-1:100
z<-runif(100)
y<-2+1.5*pmax(x-35,0)-1.5*pmax(x-70,0)+10*pmax(z-.5,0)+rnorm(100,0,2)
out.lm<-lm(y~x)
o<-segmented(out.lm,seg.Z=~x+z,psi=list(x=c(30,60),z=.4))
confint(o)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
