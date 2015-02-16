\name{seg.lm.fit}
\alias{seg.lm.fit}
\alias{seg.glm.fit}
\alias{seg.def.fit}
\alias{seg.Ar.fit}
\alias{seg.lm.fit.boot}
\alias{seg.glm.fit.boot}
\alias{seg.def.fit.boot}
\alias{seg.Ar.fit.boot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitter Functions for Segmented Linear Models
}
\description{
\code{seg.lm.fit} is called by \code{segmented.lm} to fit segmented linear 
(gaussian) models.  Likewise, \code{seg.glm.fit} is called by \code{segmented.glm} to fit  
generalized segmented linear models, and \code{seg.def.fit} is called by \code{segmented.default} to fit  
segmented relationships in general regression models (e.g., quantile regression and Cox regression). \code{seg.lm.fit.boot}, 
\code{seg.glm.fit.boot}, and \code{seg.def.fit.boot} are employed to perform bootstrap restart. 
These functions should usually not be used directly by the user.
}
\usage{
seg.lm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=FALSE)

seg.lm.fit.boot(y, XREG, Z, PSI, w, offs, opz, n.boot=10,
    size.boot=NULL, jt=FALSE, nonParam=TRUE, random=FALSE)

seg.glm.fit(y, XREG, Z, PSI, w, offs, opz, return.all.sol=FALSE)

seg.glm.fit.boot(y, XREG, Z, PSI, w, offs, opz, n.boot=10,
    size.boot=NULL, jt=FALSE, nonParam=TRUE, random=FALSE)

seg.def.fit(obj, Z, PSI, mfExt, opz, return.all.sol=FALSE)

seg.def.fit.boot(obj, Z, PSI, mfExt, opz, n.boot=10, size.boot=NULL, 
    jt=FALSE, nonParam=TRUE, random=FALSE)

seg.Ar.fit(obj, XREG, Z, PSI, opz, return.all.sol=FALSE)

seg.Ar.fit.boot(obj, XREG, Z, PSI, opz, n.boot=10, size.boot=NULL, jt=FALSE,
    nonParam=TRUE, random=FALSE)

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
  \item{offs}{
  possibe offset vector.
}
  \item{opz}{
  a list including information useful for model fitting.
}
  \item{n.boot}{
  the number of bootstrap samples employed in the bootstrap restart algorithm.
}
  \item{size.boot}{
  the size of the bootstrap resamples. If \code{NULL} (default), it is taken equal to the sample size.
  values smaller than the sample size are expected to increase perturbation in the bootstrap resamples.
}
  \item{jt}{
  logical. If \code{TRUE} the values of the segmented variable(s) are jittered before fitting the model to the
  bootstrap resamples.
}
  \item{nonParam}{
  if \code{TRUE} nonparametric bootstrap (i.e. case-resampling) is used, otherwise residual-based.

}
  \item{random}{
  if \code{TRUE}, when the algorithm fails to obtain a solution, random values are used as candidate values.
}
  \item{return.all.sol}{
  if \code{TRUE}, when the algorithm fails to obtain a solution, the values visited by the algorithm
  with corresponding deviances are returned.
}
  \item{obj}{
  the starting regression model where the segmented relationships have to be added.
}
  \item{mfExt}{
  the model frame.
}

}
\details{
The functions call iteratively \code{lm.wfit} (or \code{glm.fit}) with proper design matrix depending on 
\code{XREG}, \code{Z} and \code{PSI}. \code{seg.lm.fit.boot} (and \code{seg.glm.fit.boot}) implements the bootstrap restarting idea discussed in
Wood (2001).
}
\value{
A list of fit information.
}
\references{ Wood, S. N. (2001) Minimizing model fitting objectives that contain spurious local minima
    by bootstrap restarting. \emph{Biometrics} \bold{57}, 240--244. }
\author{ Vito Muggeo }
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
