\name{plot.segmented}
\alias{plot.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot method for segmented objects  }
\description{
  Takes a fitted \code{segmented} object returned by \code{segmented()} and plots (or adds)
  the fitted broken-line for the selected segmented term.
}
\usage{
\method{plot}{segmented}(x, term, add=FALSE, res=FALSE, se=FALSE, show.gap=TRUE, linkinv = FALSE,
      res.col=1,  rev.sgn=FALSE, const=0, ...)
      }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a fitted \code{segmented} object. }
  \item{term}{ the segmented variable having the piece-wise relationship to be plotted.
      If there is a single segmented variable, \code{term} can be omitted.}
  \item{add}{ when \code{TRUE} the fitted lines are added to the current device.}
  \item{res}{ when \code{TRUE} the fitted lines are plotted along with corresponding partial residuals.
      See Details.}
  \item{se}{ when \code{TRUE} pointwise confidence intervals are drawn. Currently unimplemented. }
  \item{show.gap}{ when \code{FALSE} the (possible) gaps between the fitted lines at the estimated breakpoints
      are hidden.}
  \item{linkinv}{ when \code{TRUE}, the fitted lines are tranformed on the inverse link (i.e. response) scale
    before plotting.}
  \item{res.col}{when \code{res=TRUE} it means the color of the points representing the partial residuals.}
  \item{rev.sgn}{ when \code{TRUE} it is assumed that current \code{term} is `minus' the actual segmented variable,
    therefore the sign is reversed before plotting. This is useful when a null-constraint has been set on the last slope.}
  \item{const}{ constant to add to each fitted segmented relationship (on the scale of the linear predictor) before plotting.}
  \item{\dots}{ other graphics parameters to pass on to plotting commands: `col', `lwd' and `lty' (that
  can be vectors, see the example below) for the fitted piecewise lines; `ylab', `xlab', `main', `sub' when a new plot
  is produced (i.e. when \code{add=FALSE}); `pch' and `cex' for the partial residuals (when \code{res=TRUE}).
  }
}
\details{
  Produces (or adds to the current device) the fitted segmented relationship between the
  response and the selected \code{term}. If the fitted model includes just a single `segmented' variable,
  \code{term} may be omitted. Due to the parameterization of the segmented terms, sometimes
  the fitted lines may not appear to join at the estimated breakpoints. If this is the case, the apparent
  `gap' would indicate some lack-of-fit. However such `gap' may be hidden by setting \code{show.gap=FALSE}:
  \code{show.gap=TRUE} is suggested only when the gap coefficients are nonsignificant, see
  \code{summary.segmented(..)$gap}.
  The partial residuals are computed as `fitted + residuals', where `fitted' are the fitted values of the
  segmented relationship. Notice that for GLMs the residuals are the response residuals if \code{linkinv=TRUE} and
  the working residuals weighted by the IWLS weights if \code{linkinv=FALSE}.
  }
\value{
None.
}
%\references{ }
\author{ Vito M.R. Muggeo }
%\note{ ~~further notes~~
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{lines.segmented}} }
\examples{
set.seed(1234)
z<-runif(100)
y<-rpois(100,exp(2+1.8*pmax(z-.6,0)))
o<-glm(y~z,family=poisson)
o.seg<-segmented(o,seg.Z=~z,psi=list(z=.5))
par(mfrow=c(1,2))
plot(o.seg)
plot(z,y)
## add the fitted lines using different colors and styles..
plot(o.seg,add=TRUE,linkinv=TRUE,lwd=2,col=2:3, lty=c(1,3))
lines(o.seg,col=2,pch=19,bottom=FALSE,lwd=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
