\name{vcov.segmented.lme}
\alias{vcov.segmented.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance-Covariance Matrix for a Fitted Segmented Mixed Model}
\description{
Returns the variance-covariance matrix of the parameters (including breakpoints) of a
     fitted segmented mixed model object.}
\usage{
\method{vcov}{segmented.lme}(object, B=0, ret.b=FALSE, ...) %is???
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted model object of class "segmented.lme", returned by \code{segmented.lme} method.}
  \item{B}{number of bootstrap replicates, if a bootstrap-based covariance matrix is requested.}
  \item{ret.b}{logical. If \code{FALSE} the full covariance matrix (for the fixed effect estimates) based on \code{B} case-resampling bootstrap samples is returned; otherwise a list with information on the bootstrap sampling distributions.}
  \item{...}{optional arguments, i.e. \code{seed} and \code{it.max.b}, used when implementing the bootstrap.}
%  \item{is}{logical. If \code{TRUE}, the \emph{asymptotic} covariance matrix based on the idea of induced smoothing is returned. If \code{is%=TRUE}, \code{var.diff=FALSE} is set. 
%  }
}
\details{
The returned covariance matrix is based on an approximation of the nonlinear segmented term. Therefore 
covariances corresponding to breakpoints are reliable only in large samples and/or clear cut segmented 
relationships. If \code{B>0} is set, case resampling bootstrap (on the outermost nesting level) is carried out. Moreover, if \code{ret.b=TRUE}, the bootstrap distributions are returned, rather than the covariance matrix.  
%If \code{is=TRUE}, the returned covariance matrix depends on the design matrix having the term \eqn{I(x>\psi)}{I(x>psi)} replaced by its smooth %counterpart.
}
\value{
The full matrix of the estimated covariances of the fixed effects estimates, including 
the breakpoint.
}
%\references{}
\author{Vito M. R. Muggeo, \email{vito.muggeo@unipa.it}}

\section{Warning }{
All the functions for segmented mixed models (*.segmented.lme) are still at an experimental stage 
} 

\seealso{\code{\link{summary.segmented.lme}}}

\examples{
##continues example from segmented.lme()
# vcov(os)
# vcov(os, B=50)
# vcov(os, B=50, ret.b=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
