\name{vcov.stepmented}
\alias{vcov.stepmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance-Covariance Matrix for a Fitted Stepmented Model
}
\description{
Returns the variance-covariance matrix of the parameters estimates (including breakpoints) of a
     fitted stepmented model object.
}
\usage{
\method{vcov}{stepmented}(object, k=NULL, zero.cor=TRUE, type=c("cdf", "none", "abs"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted model object of class "stepmented", returned by any \code{stepmented} method
}
  \item{k}{
The power of \code{n} for the smooth approximation. Simulation evidence suggests \code{k} in \eqn{[-1, -1/2]}; with \eqn{k=-1/2} providing somewhat 'conservative' standard errors especially at small sample sizes. In general, the larger \eqn{k}, the smaller \eqn{n^{-k}}, and the smaller the jumpoint standard error.  
}
\item{zero.cor}{
If \code{TRUE}, the covariances between the jumpoints and the remaining linear coefficients are set to zero (as theory states).
}
\item{type}{How the covariance matrix should be computed. If \code{"none"}, the usual asymptotic covariance matrix for the linear coefficients only (under homoskedasticity and assuming known the jumpoints) is returned; if \code{"cdf"}, the standard normal cdf is used to approximate the indicator function (see details); \code{"abs"} is yet another approximation (currently unimplemented). 
}
\item{\dots}{
additional arguments. 
}
}
\details{
The full covariance matrix is based on the smooth approximation 
\deqn{I(x>\psi)\approx \Phi((x-\psi)/n^{k})} 
via the sandwich formula using the empirical information matrix and assuming \eqn{x \in [0,1]}. \eqn{\Phi(\cdot)} is the standard Normal cdf, and \eqn{k} is the argument \code{k}. When \code{k=NULL} (default), it is computed via
\deqn{k=-(0.6 + 0.5 \ \log(snr)/\sqrt snr - (|\hat\psi-0.5|/n)^{1/2})}
 where \eqn{snr} is the signal-to-noise ratio corresponding to the estimated changepoint \eqn{\hat\psi} (in the range (0,1)). The above formula comes from extensive simulation studies under different scenarios: Seo and Linton (2007) discuss using the normal cdf to smooth out the indicator function by suggesting \eqn{\log(n)/n^{1/2}} as bandwidth; we found such suggestion does not perform well in practice. 
}
\value{
The full matrix of the estimated covariances between the parameter estimates, including 
the breakpoints.
}
\references{
Seo MH, Linton O (2007) A smoothed least squares estimator for threshold regression models, J of Econometrics, 141: 704-735
}
\author{
Vito Muggeo
}
\note{
If the fit \code{object} has been called by \code{stepmented(.., var.psi=TRUE)}, then \code{vcov.stepmented} will return \code{object$vcov}, unless the power \code{k} differs from \code{-2/3}.  
}

\section{Warning }{The function, including the value of \eqn{k}, must be considered at preliminary stage. Currently the value of \eqn{k} appears to  overestimate slightly the true \eqn{\hat\psi} variability.} 

\seealso{
\code{\link{stepmented}}
}
\examples{
##see ?stepmented
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
