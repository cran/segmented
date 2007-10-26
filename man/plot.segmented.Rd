\name{plot.segmented}
\alias{plot.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot method for segmented objects  }
\description{
  Takes a fitted \code{segmented} object returned by \code{segmented()} and plots (or add)
  the fitted broken-line for the selected segmented term.
}
\usage{
\method{plot}{segmented}(x, term = NULL, se = FALSE, const = coef(x)["(Intercept)"],
  add = FALSE, linkinv = FALSE, h = 300, show.gap=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a fitted \code{segmented} object }
  \item{term}{ the segmented variable whose piece-wise relationship has to be plotted }
  \item{se}{ when \code{TRUE} pointwise confidence intervals are drawn. Currently unimplemented. }
  \item{const}{ constant to add to each fitted segmented relationship (on the scale of the linear predictor) before plotting }
  \item{add}{ when \code{TRUE} the fitted lines are added to the current device }
  \item{linkinv}{ when \code{TRUE}, the fitted lines are (possibly) tranformed on the inverse link scale
    before plotting}
  \item{h}{ graphic parameter to regulate the number of points where the fitted lines have to be computed.
    See Details  }
  \item{show.gap}{ should the gap between the fitted lines at the estimated breakpoints to be shown?}
  \item{\dots}{ other graphics parameters to pass on to plotting commands }
}
\details{
  Produces (or adds to the current device) the fitted segmented relationship between the
  response and the selected \code{term}. If the fitted model includes just a single `segmented' variable,
  \code{term} may be omitted. Due to the parameterization of the segmented terms, sometimes
  the fitted lines may not appear to join at the estimated breakpoints. If this is the case, try to
  increase the value of \code{h}; otherwise the apparent `gap' would indicate some lack-of-fit.
  However such `gap' may be hidden by setting \code{show.gap=FALSE}: in this case the new fitted
  values are re-computed by means of a linear model fit. Note that, as hidding could lead to a wrong
  impression of the fit, \code{show.gap=TRUE} is suggested only when the gap coefficients are nonsignificant.
  }
\value{
None.
}
\references{ }
\author{ Vito M.R. Muggeo }
%\note{ ~~further notes~~
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{point.psi}} }
\examples{
set.seed(1234)
z<-runif(100)
y<-rpois(100,exp(2+1.8*pmax(z-.6,0)))
o<-glm(y~z,family=poisson)
o.seg<-segmented(o,seg.Z=~z,psi=list(z=.5))
par(mfrow=c(1,2))
plot(o.seg)
plot(z,y)
plot(o.seg,add=TRUE,linkinv=TRUE,lwd=2,col=2)
point.psi(o.seg,"z",col=2,pch=19,bottom=FALSE,lwd=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
