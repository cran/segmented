\name{plot.stepmented}
\alias{plot.stepmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for stepmented objects
}
\description{
  Takes a fitted \code{stepmented} object returned by \code{stepmented()} and plots (or adds)
  the fitted piecewise constant lines for the selected stepmented term.
}
\usage{
\method{plot}{stepmented}(x, term, add = FALSE, res = TRUE, conf.level=0, interc = TRUE, 
    add.fx = FALSE, psi.lines = TRUE, link=FALSE, surf=FALSE, zero.cor=TRUE, 
                            heurs=TRUE, shade=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a fitted \code{stepmented} object.
}
  \item{term}{
the stepmented variable having the piece-wise constant relationship to be plotted.
If there is a single stepmented variable in the fitted model \code{x}, 
\code{term} can be omitted.
}
  \item{add}{
when \code{TRUE} the fitted lines are added to the current device.
}
  \item{res}{
when \code{TRUE} the fitted lines are plotted along with corresponding partial residuals.
}
  \item{conf.level}{
the confidence level for the pointwise confidence intervals for the expected values.
}

  \item{interc}{
if \code{TRUE} the computed components include the model intercept (if it exists).
}
  \item{add.fx}{
logical. If TRUE and the object fit also includes an additional term for the same stepmented variable, the plot also portrays such `additional' term. 
}
  \item{psi.lines}{
if \code{TRUE} vertical lines corresponding to the estimated changepoints are also drawn
}
  \item{link}{
if \code{FALSE} the fitted lines (and possibily the residuals) are reported on the response scale. Ignored if the fit object \code{x} is not a glm-like fit.
}
  \item{surf}{
if the object fit \code{x} includes 2 stepmented covariates (x1 and x2, say) with relevant estimated breakpoints, \code{surf=TRUE} will draw on the plane x1-x2 the areas splitted according to the estimated breakpoints with corresponding estimated means superimposed.
}
  \item{zero.cor}{
see \code{zero.cor} in \code{\link{vcov.stepmented}}; effective only if \code{conf.level>0}.
}
  \item{heurs}{
logical; if \code{TRUE}, heuristic (usually somewhat conservative) confidence intervals are computed and plotted; effective only if \code{conf.level>0}. 
}
  \item{shade}{
if \code{TRUE} the pointwise confidence intervals are portrayed via shaded area; effective only if \code{conf.level>0}.
}

  \item{\dots}{ other graphics parameters to pass to plotting commands: `col', `lwd' and `lty' (that
  can be vectors and are recycled if necessary, see the example below) for the fitted piecewise constant lines; `ylab', `xlab', `main', `sub', `cex.axis', `cex.lab', `xlim' and `ylim' when a new plot is produced (i.e. when \code{add=FALSE}); `pch' and `cex' for the partial residuals (when \code{res=TRUE}, \code{res.col} is for the color).% \code{col.shade} for the shaded regions (provided that \code{shade=TRUE} and \code{conf.level>0}).
  }
}
\details{
  Produces (or adds to the current device) the fitted step-function like relationship between the
  response and the selected \code{term}. If the fitted model includes just a single `stepmented' variable,
  \code{term} may be omitted. If \code{surf=TRUE}, and \code{res=TRUE} the point widths are proportional to the partial residual values. 
}
\value{
None.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M. R. Muggeo
}
\note{
Implementation of confidence intervals for the conditional means in stepmented regression is under development; \code{conf.level>0} should be used with care, especially with multiple jumpoints.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See Also as \code{\link{stepmented}}
}
\examples{

#Following code in stepmented..
\dontrun{

par(mfrow=c(1,3))

plot(os,"x")
plot(os,"z")
plot(os,"z", add.fx=TRUE, psi.lines=FALSE ) 
lines(os, "z")

#display the 'surface'
par(mfrow=c(1,2))
plot(os, surf=TRUE, col=1, res.col=2)
plot(os, surf=TRUE, lty=2)
}

}
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ hplot }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
