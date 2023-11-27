\name{plot.segmented}
\alias{plot.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot method for segmented objects  }
\description{
  Takes a fitted \code{segmented} object returned by \code{segmented()} and plots (or adds)
  the fitted broken-line relationship for the selected segmented term.
}
\usage{
\method{plot}{segmented}(x, term, add=FALSE, res=FALSE, conf.level=0, interc=TRUE, 
    link=TRUE, res.col=grey(.15, alpha = .4), rev.sgn=FALSE, const=NULL, shade=FALSE, 
    rug=!add, dens.rug=FALSE, dens.col = grey(0.8), transf=I, 
    isV=FALSE, is=FALSE, var.diff=FALSE, p.df="p", .vcov=NULL, .coef=NULL, 
    prev.trend=FALSE, smoos=NULL, hide.zeros=FALSE, ...)
      }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a fitted \code{segmented} object. }
  \item{term}{ the segmented variable having the piece-wise relationship to be plotted.
      If there is a single segmented variable in the fitted model \code{x}, \code{term} can be omitted.}
  \item{add}{ when \code{TRUE} the fitted lines are added to the current device.}
  \item{res}{ when \code{TRUE} the fitted lines are plotted along with corresponding partial residuals.
      See Details.}
  \item{conf.level}{ If greater than zero, it means the confidence level at which the pointwise confidence itervals have to 
  be plotted.}
  \item{interc}{ If \code{TRUE} the computed segmented components include the model intercept (if it exists).}
  \item{link}{ when \code{TRUE} (default), the fitted lines are plotted on the link scale, otherwise they are
    tranformed on the response scale before plotting. Ignored for linear segmented fits. }
  \item{res.col}{when \code{res=TRUE} it means the color of the points representing the partial residuals.}
  \item{rev.sgn}{ when \code{TRUE} it is assumed that current \code{term} is `minus' the actual segmented variable,
    therefore the sign is reversed before plotting. This is useful when a null-constraint has been set on the last slope.}
  \item{const}{ constant to add to each fitted segmented relationship (on the scale of the linear predictor) before plotting.
  If the fit includes a segmented interaction term (obtained via \code{seg(..,by)} in the formula), the group-specific intercept is included. 
  }
  \item{shade}{if \code{TRUE} and \code{conf.level>0} it produces shaded regions (in grey color) for the pointwise confidence 
  intervals embracing the fitted segmented line. }
  \item{rug}{when \code{TRUE} the covariate values are displayed as a rug plot
 at the foot of the plot. Default is to \code{!add}.}
  \item{dens.rug}{when \code{TRUE} then smooth covariate distribution is plotted on the x-axis.}
  \item{dens.col}{if \code{dens.rug=TRUE}, it means the colour to be used to plot the density.}
%  \item{show.gap}{ when \code{FALSE} the (possible) gaps between the fitted lines at the estimated breakpoints
%      are hidden. When bootstrap restarting has been employed (default in \code{segmented}), \code{show.gap} is meaningless 
%      as the gap coefficients are always set to zero in the fitted model.}
  \item{transf}{ A possible function to convert the fitted values before plotting. It is only effective 
    if \code{res=FALSE}. If \code{res=TRUE} any transformation is ignored.}
  \item{isV}{logical value (to be passed to \code{\link{broken.line}}). Ignored if \code{conf.level=0}}
  \item{is}{logical value (to be passed to \code{\link{broken.line}}) indicating if the covariance matrix based on the induced smoothing should be used. Ignored if \code{conf.level=0}}
  \item{var.diff}{logical value to be passed to \code{\link{summary.segmented}} to compute dthe standard errors of fitted values (if \code{conf.level>0}).}
  \item{p.df}{ degrees of freedom when \code{var.diff=TRUE}, see \code{\link{summary.segmented}}}
  \item{.vcov}{ The \emph{full} covariance matrix of estimates to be used when \code{conf.level>0}. If unspecified (i.e. \code{NULL}), the covariance matrix is computed internally by the function \code{vcov.segmented}.}
  \item{.coef}{ The regression parameter estimates. If unspecified (i.e. \code{NULL}), it is computed internally by \code{coef()}.}
  \item{prev.trend}{logical. If \code{TRUE} dashed lines corresponding to the `previous' trends (i.e. the trends if the breakpoints would not have occurred) are also drawn.}
    \item{smoos}{
logical, indicating if the residuals (provided that \code{res=TRUE}) will be drawn using a \emph{smoothed} scatterplot. If \code{NULL} (default) the smoothed scatterplot will be employed when the number of observation is larger than 10000.
}
    \item{hide.zeros}{
logical, indicating if the residuals (provided that \code{res=TRUE}) corresponding to the covariate zero values should be deleted. Useful when the fit includes an interaction term in the formula, such as \code{seg(.., by=..)}, and the zeros in covariates indicate units in other groups.
}
  \item{\dots}{ other graphics parameters to pass to plotting commands: `col', `lwd' and `lty' (that
  can be vectors and are recycled if necessary, see the example below) for the fitted piecewise lines; `ylab', `xlab', `main', `sub', `cex.axis', `cex.lab', `xlim' and 
  `ylim' when a new plot is produced (i.e. when \code{add=FALSE}); `pch' and `cex' for the partial residuals (when \code{res=TRUE}, \code{res.col} is for the color); \code{col.shade} for the shaded regions (provided that \code{shade=TRUE} and \code{conf.level>0}).
  }
}
\details{
  Produces (or adds to the current device) the fitted segmented relationship between the
  response and the selected \code{term}. If the fitted model includes just a single `segmented' variable,
  \code{term} may be omitted. 
  %Due to the parameterization of the segmented terms, sometimes
  %the fitted lines may not appear to join at the estimated breakpoints. If this is the case, the apparent
  %`gap' would indicate some lack-of-fit. However, since version 0.2-9.0, the gap coefficients are set to zero by default 
  %(see argument \code{gap} in in \code{\link{seg.control}}). 
  The partial residuals are computed as `fitted + residuals', where `fitted' are the fitted values of the
  segmented relationship relevant to the covariate specified in \code{term}. 
  Notice that for GLMs the residuals are the response residuals if \code{link=FALSE} and
  the working residuals if \code{link=TRUE}. %weighted by the IWLS weights  [fino alla versione 0.5-2.0 i workRes were weighted by the IWLS weights]
  }
\value{
None.
}
%\references{ }
\author{ Vito M. R. Muggeo }
\note{ For models with offset, partial residuals on the response scale are not defined. Thus \code{plot.segmented} does not work 
when \code{link=FALSE}, \code{res=TRUE}, and the fitted model includes an offset.}
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{segmented}} to fit the model, \code{\link{lines.segmented}} to add the estimated breakpoints on the current plot. 
\code{\link{points.segmented}} to add the joinpoints of the segmented relationship.
\code{\link{predict.segmented}} to compute standard errors and confidence intervals for predictions from a 
"segmented" fit. }
\examples{
set.seed(1234)
z<-runif(100)
y<-rpois(100,exp(2+1.8*pmax(z-.6,0)))
o<-glm(y~z,family=poisson)
o.seg<-segmented(o) #single segmented covariate and one breakpoint: 'seg.Z' and 'psi' can be omitted
par(mfrow=c(1,2))
plot(o.seg, conf.level=0.95, shade=TRUE)
points(o.seg, link=TRUE, col=2)
## new plot
plot(z,y)
## add the fitted lines using different colors and styles..
plot(o.seg,add=TRUE,link=FALSE,lwd=2,col=2:3, lty=c(1,3))
lines(o.seg,col=2,pch=19,bottom=FALSE,lwd=2) #for the CI for the breakpoint
points(o.seg,col=4, link=FALSE)
## using the options 'is', 'isV', 'shade' and 'col.shade'.
par(mfrow=c(1,2))
plot(o.seg, conf.level=.9, is=TRUE, isV=TRUE, col=1, shade = TRUE, col.shade=2)
plot(o.seg, conf.level=.9, is=TRUE, isV=FALSE, col=2, shade = TRUE, res=TRUE, res.col=4, pch=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ hplot }