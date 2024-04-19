\name{confint.stepmented}
\alias{confint.stepmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence intervals for jumpoints in stepmented regression}
\description{
  Computes confidence intervals for the changepoints (or jumpoints) in a fitted `stepmented' model.
}
\usage{
\method{confint}{stepmented}(object, parm, level=0.95, method=c("delta", "score", "gradient"), %var.diff=FALSE, 
    round=TRUE, cheb=FALSE, digits=max(4, getOption("digits") - 1), 
    .coef=NULL, .vcov=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted \code{stepmented} object. }
  \item{parm}{the stepmented variable of interest. If missing the first stepmented variable in \code{object} is considered. }
  \item{level}{the confidence level required, default to 0.95.}
  \item{method}{which confidence interval should be computed. One of \code{"delta"}, \code{"score"}, or \code{"gradient"}. Can be abbreviated. Currently only \code{"delta"} allowed.}
  \item{round}{logical. Should the values (estimates and lower/upper limits) rounded to the smallest observed value?}
  \item{cheb}{logical. If \code{TRUE}, the confidence limits are computed using the Chebyshev inequality which yields conservative confidence intervals but it is 'robust' to the non-normality of the changepoint sampling distribution. 
  }
  %\item{var.diff}{logical. If \code{method="delta"},  and there is a single segmented variable, \code{var.diff=TRUE} leads     to standard errors based on sandwich-type formula of the covariance matrix. See Details in \code{\link{summary.segmented}}.}
  %\item{is}{logical. If \code{method="delta"}, \code{is=TRUE} means that the full covariance matrix is computed via 
  %\code{vcov(.., is=TRUE)}}
  \item{digits}{controls the number of digits to print when returning the output. }
  \item{.coef}{ The regression parameter estimates. If unspecified (i.e. \code{NULL}), it is computed internally by \code{coef(object)}.}
  \item{.vcov}{ The \emph{full} covariance matrix of estimates. If unspecified (i.e. \code{NULL}), the covariance matrix is computed internally by \code{vcov(object)}.}
  \item{\dots}{additional arguments passed to \code{vcov.stepmented}, namely \code{k}. 
  %referring to Score-based confidence intervals, such as \code{"h"}, \code{"d.h"}, \code{"bw"}, \code{"msgWarn"}, and \code{"n%.values"} specifying the number of points used to profile the Score (or Gradient) statistic.
  }
}
\details{
  \code{confint.stepmented} computes confidence limits for the changepoints. Currently the only option is \code{'delta'}, i.e. to compute the approximate covariance matrix via a smoothing approximation (see \code{\link{vcov.stepmented}}) and to build the limits using the standard Normal quantiles. Note that, the limits are rounded to the lowest observed value, thus the resulting confidence interval might not be symmetric if the stepmented covariate has not equispaced values. 
  }
\value{
    A matrix including point estimate and confidence limits of the breakpoint(s) for the
    stepmented variable possibly specified in \code{parm}.
}
%\references{ 
%Muggeo, V.M.R. (2017) Interval estimation for the breakpoint in segmented regression: a smoothed score-based approach. 
%\emph{Australian & New Zealand Journal of Statistics} \bold{59}, 311--322.
%}
\author{ Vito M.R. Muggeo }
\note{ 
Currently only method='delta' is allowed. 
%\code{method="score"} or \code{method="gradient"} only works for segmented \emph{linear} model. For segmented \emph{generalized linear} model, 
 % currently only \code{method="delta"} is available. 
}
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{stepmented}} and \code{\link{lines.segmented}} to plot the estimated breakpoints with corresponding 
confidence intervals. }
\examples{
set.seed(10)
x<-1:100
z<-runif(100)
y<-2+2.5*(x>45)-1.5*(x>70)+z+rnorm(100)
o<-stepmented(y, npsi=2)

confint(o) #round=TRUE is default
confint(o, round=FALSE)  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
