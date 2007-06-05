\name{segmented.default}
\alias{segmented.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Default method for the generic segmented }
\description{
  \code{segmented} is a generic function, and \code{segmented.default} its default method.
  }
\usage{
segmented.default(obj, seg.Z, psi, control = seg.control(), model.frame = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ see \code{\link{segmented}} }
  \item{seg.Z}{ see \code{\link{segmented}} }
  \item{psi}{ see \code{\link{segmented}} }
  \item{control}{ see \code{\link{segmented}} }
  \item{model.frame}{ see \code{\link{segmented}} }
  \item{\dots}{ see \code{\link{segmented}} }
}
\details{
  Actually \code{segmented.default} makes nothing! Use the specific methods.
}
\seealso{ \code{\link{segmented}} }
\examples{
##Does not work!
# segmented.default(obj.glm,..)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
