\name{points.segmented}
\alias{points.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Points method for segmented objects
}
\description{
  Takes a fitted \code{segmented} object returned by \code{segmented()} and adds
  on the current plot the joinpoints of the fitted broken-line relationships.
}
\usage{
  \method{points}{segmented}(x, term, interc = TRUE, link = TRUE, rev.sgn=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{segmented}.
}
  \item{term}{
the segmented variable of interest. It may be unspecified when there is a single segmented variable.
}
  \item{interc}{
If \code{TRUE} the computed joinpoints include the model intercept (if it exists).
}
  \item{link}{
when \code{TRUE} (default), the fitted joinpoints are plotted on the link scale
}
  \item{rev.sgn}{
when \code{TRUE}, the fitted joinpoints are plotted on the `minus' scale
of the current \code{term} variable. This is useful when a null-constraint 
has been set on the last slope.
}
  \item{\dots}{
other graphics parameters to pass on to \code{points()} function.
}
}
\details{
We call 'joinpoint' the plane point having as coordinates the breakpoint (on the x scale) and the fitted value of
the segmented relationship at that breakpoint (on the y scale). \code{points.segmented()} simply adds the fitted 
joinpoints on the current plot. This could be useful to emphasize the changes of the piecewise linear relationship. 
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.segmented}} to plot the fitted segmented lines. }

\examples{
\dontrun{
#continues from ?plot.segmented
points(o.seg,col=2)
}
}
\keyword{ nonlinear }
\keyword{ regression }% __ONLY ONE__ keyword per line
