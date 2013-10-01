\name{draw.history}
\alias{draw.history}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ History for the breakpoint estimates }
\description{
  Displays breakpoint iteration values for segmented fits.
}
\usage{
draw.history(obj, term, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ a segmented fit returned by any "segmented" method. }
  \item{term}{ a character to mean the `segmented' variable whose breakpoint values throughout iterations have to be displayed. }
  \item{\dots}{ graphic parameters to be passed to \code{matplot()}. }
}
\details{
  For a given \code{term} in a segmented fit, \code{draw.history()} displays the different breakpoint
  values obtained during the estimating process, since the starting values up to the final ones. When 
  bootstrap restarting is employed, \code{draw.history()} produces two plots, the values of objective function 
  and the number of distinct solutions against the bootstrap replicates.
}
\value{ None. }
%\references{  }
\author{ Vito M.R. Muggeo }
%\note{ ~~further notes~~
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(stagnant)
os<-segmented(lm(y~x,data=stagnant),seg.Z=~x,psi=-.8)
draw.history(os) #diagnostics with boot restarting

os<-segmented(lm(y~x,data=stagnant),seg.Z=~x,psi=-.8, control=seg.control(n.boot=0))
draw.history(os) #diagnostics without boot restarting

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }

