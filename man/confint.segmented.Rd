\name{confint.segmented}
\alias{confint.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence intervals for breakpoints}
\description{
  Computes confidence intervals for the breakpoints in a fitted `segmented' model.
}
\usage{
\method{confint}{segmented}(object, parm, level=0.95, method=c("delta", "score", "gradient"), 
    rev.sgn=FALSE, var.diff=FALSE, is=FALSE, digits=max(4, getOption("digits") - 1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted \code{segmented} object. }
  \item{parm}{the segmented variable of interest. If missing the first segmented variable in \code{object} is considered. }
  \item{level}{the confidence level required, default to 0.95.}
  \item{method}{which confidence interval should be computed. One of \code{"delta"}, \code{"score"}, or \code{"gradient"}. Can be abbreviated.}
  \item{rev.sgn}{vector of logicals. The length should be equal to the length of \code{parm}; recycled otherwise.
   when \code{TRUE} it is assumed that the current \code{parm} is `minus' the actual segmented variable,
    therefore the sign is reversed before printing. This is useful when a null-constraint has been set on the last slope.}
  \item{var.diff}{logical. If \code{method="delta"},  and there is a single segmented variable, \code{var.diff=TRUE} leads     to standard errors based on sandwich-type formula of the covariance matrix. See Details in \code{\link{summary.segmented}}.}
  \item{is}{logical. If \code{method="delta"}, \code{is=TRUE} means that the full covariance matrix is computed via 
  \code{vcov(.., is=TRUE)}}
  \item{digits}{controls the number of digits to print when returning the output. }
  \item{\dots}{additional parameters referring to Score-based confidence intervals, such as \code{"h"}, \code{"d.h"}, \code{"bw"}, \code{"msgWarn"}, and \code{"n.values"} 
  specifying the number of points used to profile the Score (or Gradient) statistic.}
}
\details{
  \code{confint.segmented} computes confidence limits for the breakpoints. Currently there are three options, see argument \code{method}. 
  \code{method="delta"} uses the standard error coming from the Delta
  method for the ratio of two random variables. This value is an approximation (slightly) better than the 
  one reported in the `psi' component of the list returned by any \code{segmented} method. The resulting
  confidence intervals are based on the asymptotic Normal distribution of the breakpoint
  estimator which is reliable just for clear-cut kink relationships. See Details in \code{\link{segmented}}. \cr
  \code{method="score"} or \code{method="gradient"} compute the 
  confidence interval via profiling the Score or the Gradient statistics smoothed out by the induced smoothing paradigm, as discussed in the reference below. 
  }
\value{
    A matrix including point estimate and confidence limits of the breakpoint(s) for the
    segmented variable possibly specified in \code{parm}.
}
\references{ 
Muggeo, V.M.R. (2017) Interval estimation for the breakpoint in segmented regression: a smoothed score-based approach. 
\emph{Australian & New Zealand Journal of Statistics} \bold{59}, 311--322.
}
\author{ Vito M.R. Muggeo }
\note{ 
Currently \code{method="score"} or \code{method="gradient"} only works for segmented \emph{linear} model. For segmented \emph{generalized linear} model, 
  currently only \code{method="delta"} is available. 
}
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{segmented}} and \code{\link{lines.segmented}} to plot the estimated breakpoints with corresponding 
confidence intervals. }
\examples{
set.seed(10)
x<-1:100
z<-runif(100)
y<-2+1.5*pmax(x-35,0)-1.5*pmax(x-70,0)+10*pmax(z-.5,0)+rnorm(100,0,2)
out.lm<-lm(y~x)
o<-segmented(out.lm,seg.Z=~x+z,psi=list(x=c(30,60),z=.4))
confint(o) #delta CI for the 1st variable 
confint(o, "x", method="score") #also method="g"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
