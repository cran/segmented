\name{seg.lm.fit}
\alias{seg.lm.fit}
\alias{seg.glm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitter Function for Segmented Linear Models
}
\description{
\code{seg.lm.fit} is called by \code{segmented.lm} to fit segmented linear 
(gaussian) models.  Likewise \code{seg.glm.fit} is called by \code{segmented.glm} to fit  
generalized segmented linear models. These functions should usually not be used directly by the user.
}
\usage{
seg.lm.fit(y, XREG, Z, PSI, w, o, opz)

seg.glm.fit(y, XREG, Z, PSI, w, o, opz)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
vector of observations of length \code{n}.
}
  \item{XREG}{
design matrix for standard linear terms.
}
  \item{Z}{
appropriate matrix including the segmented variables whose breakpoints have to be estimated.
}
  \item{PSI}{
  appropriate matrix including the starting values of the breakpoints to be estimated.
}
  \item{w}{
  possibe weights vector.
}
  \item{o}{
  possibe offset vector.
}
  \item{opz}{
  a list including information useful for model fitting.
}
}
\details{
The functions call iteratively \code{lm.wfit} (or \code{glm.fit}) with proper design matrix depending on 
\code{XREG}, \code{Z} and \code{PSI}.
}
\value{
A list of fit information.
}
\note{
These functions should usually not be used directly by the user.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{segmented.lm}}, \code{\link{segmented.glm}}
}
\examples{
##See ?segmented
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{nonlinear }
