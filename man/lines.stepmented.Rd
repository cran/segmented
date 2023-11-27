\name{lines.stepmented}
\alias{lines.stepmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bars for interval estimate of the breakpoints }
\description{
  Draws bars relevant to breakpoint estimates (point estimate and confidence limits)
  on the current device
}
\usage{
\method{lines}{stepmented}(x, term, bottom = TRUE, shift=FALSE, conf.level = 0.95, k = 50, 
    pch = 18, .vcov=NULL, .coef=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{stepmented}. }
  \item{term}{ the stepmented variable of the breakpoints being drawn. It may be unspecified when there is a single stepmented variable.}
  \item{bottom}{ logical, indicating if the bars should be plotted at the bottom (\code{TRUE}) or at the top (\code{FALSE}).}
  \item{shift}{ logical, indicating if the bars should be `shifted' on the y-axis before plotting. Useful for multiple breakpoints
        with overlapped confidence intervals.}
  \item{conf.level}{ the confidence level of the confidence intervals for the breakpoints. }
  \item{k}{ a positive integer regulating the vertical position of the drawn bars. See Details. }
  \item{pch}{ either an integer specifying a symbol or a single character to be used
        in plotting the point estimates of the breakpoints. See \code{\link{points}}. }
    \item{.vcov}{ The \emph{full} covariance matrix of estimates. If unspecified (i.e. \code{NULL}), the covariance matrix is computed internally by \code{vcov(x)}.}
  \item{.coef}{ The regression parameter estimates. If unspecified (i.e. \code{NULL}), it is computed internally by \code{coef(x)}.}
  \item{\dots}{ further arguments passed to \code{\link{arrows}}, for instance `col' that 
  can be a vector. }
}
\details{
  \code{lines.stepmented} simply draws on the current device the point estimates and relevant
  confidence limits of the estimated breakpoints from a "stepmented" object. The y coordinates
  where the bars are drawn is computed as \code{usr[3]+h} if \code{bottom=TRUE} or
  \code{usr[4]-h} when \code{bottom=FALSE}, where \code{h=(usr[4]-usr[3])/abs(k)} and
  \code{usr} are the extremes of the user coordinates of the plotting region.
  Therefore for larger values of \code{k} the bars are plotted on the edges.
  
}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{plot.stepmented}} to plot the fitted stepmented lines}

%\code{\link{points.segmented}} to add the fitted joinpoints. }
\examples{
## See ?plot.stepmented

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
