\name{lines.segmented}
\alias{lines.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bars for interval estimate of the breakpoints }
\description{
  Draws bars relevant to breakpoint estimates (point estimate and confidence limits)
  on the current device
}
\usage{
\method{lines}{segmented}(x, term, bottom = TRUE, shift=TRUE, conf.level = 0.95, k = 50, 
    pch = 18, rev.sgn = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{segmented} }
  \item{term}{ the segmented variable of the breakpoints being drawn. It may be unspecified when there is a single segmented variable}
  \item{bottom}{ logical, indicating if the bars should be plotted at the bottom (TRUE) or at the top (FALSE)}
  \item{shift}{ logical, indicating if the bars should be `shifted' on the y-axis before plotting. Useful for multiple breakpoints
        with overlapped confidence intervals}
  \item{conf.level}{ the confidence level of the confidence intervals for the breakpoints }
  \item{k}{ a positive integer regulating the vertical position of the drawn bars. See Details }
  \item{pch}{ either an integer specifying a symbol or a single character to be used
        in plotting the point estimates of the breakpoints. See \code{\link{points}} }
  \item{rev.sgn}{ should the signs of the breakpoint estimates be changed before plotting? see Details }
  \item{\dots}{ further arguments passed to \code{\link{segments}} }
}
\details{
  \code{lines.segmented} simply draws on the current device the point estimates and relevant
  confidence limits of the estimated breakpoints from a "segmented" object. The y coordinate
  where the bars are drawn is computed as \code{usr[3]+h} if \code{bottom=TRUE} or
  \code{usr[4]-h} when \code{bottom=FALSE}, where \code{h=(usr[4]-usr[3])/abs(k)} and
  \code{usr} are the extremes of the user coordinates of the plotting region.
  Therefore for larger values of \code{k} the bars are plotted on the edges.
  The argument \code{rev.sgn} allows to change the sign of the breakpoints before plotting. This may
  be useful when a null-right-slope constraint is set.
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
\seealso{ \code{\link{plot.segmented}} }
\examples{
## See ?plot.segmented

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
