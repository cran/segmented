\name{vcov.segmented}
\alias{vcov.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance-Covariance Matrix for a Fitted Segmented Model}
\description{
Returns the variance-covariance matrix of the parameters (including breakpoints) of a
     fitted segmented model object.}
\usage{
\method{vcov}{segmented}(object, var.diff = FALSE, is = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted model object of class "segmented", returned by any \code{segmented} method.}
  \item{var.diff}{logical. If \code{var.diff=TRUE} and there is a single segmented variable, the covariance 
  matrix is computed using a sandwich-type formula. See Details in \code{\link{summary.segmented}}.}
  \item{is}{logical. If \code{TRUE}, the \emph{asymptotic} covariance matrix based on the idea of induced smoothing is returned. If \code{is=TRUE}, \code{var.diff=FALSE} is set. }
  \item{\dots}{additional arguments. }
}
\details{
The returned covariance matrix is based on an approximation of the nonlinear segmented term. Therefore 
covariances corresponding to breakpoints are reliable only in large samples and/or clear cut segmented 
relationships. If \code{is=TRUE}, the returned covariance matrix depends on the design matrix having the term \eqn{I(x>\psi)}{I(x>psi)} replaced by its smooth counterpart.
}
\value{
The full matrix of the estimated covariances between the parameter estimates, including 
the breakpoints.
}
%\references{}
\author{Vito M. R. Muggeo, \email{vito.muggeo@unipa.it}}
\note{\code{var.diff=TRUE} works when there is a single segmented variable.}
\seealso{\code{\link{summary.segmented}}}
\examples{
##continues example from summary.segmented()
# vcov(oseg)
# vcov(oseg, var.diff=TRUE)
# vcov(oseg, is=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
