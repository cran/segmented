\name{step.lm.fit}
\alias{step.lm.fit}
\alias{step.glm.fit}
\alias{step.num.fit}
\alias{step.ts.fit}
\alias{step.lm.fit.boot}
\alias{step.glm.fit.boot}
\alias{step.num.fit.boot}
\alias{step.ts.fit.boot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitter Functions for stepmented Linear Models
}
\description{
\code{step.lm.fit} is called by \code{stepmented.lm} to fit stepmented linear 
(gaussian) models.  Likewise, \code{step.glm.fit} is called by \code{stepmented.glm} to fit  
generalized stepmented linear models.
%, and \code{step.def.fit} is called by \code{stepmented.default} to fit  
%stepmented relationships in general regression models (e.g., quantile regression and Cox regression). 
The \code{step.*.fit.boot} functions are employed to perform bootstrap restarting. 
These functions should usually not be used directly by the user.
}
\usage{
step.lm.fit(y, x.lin, Xtrue, PSI, ww, offs, opz, return.all.sol=FALSE)  

step.lm.fit.boot(y, XREG, Z, PSI, w, offs, opz, n.boot=10, size.boot=NULL, 
        jt=FALSE, nonParam=TRUE, random=FALSE, break.boot=n.boot)                          

step.glm.fit(y, x.lin, Xtrue, PSI, ww, offs, opz, return.all.sol=FALSE)

step.glm.fit.boot(y, XREG, Z, PSI, w, offs, opz, n.boot=10, size.boot=NULL, 
        jt=FALSE, nonParam=TRUE, random=FALSE, break.boot=n.boot)
%step.def.fit(obj, Z, PSI, mfExt, opz, return.all.sol=FALSE)
%step.def.fit.boot(obj, Z, PSI, mfExt, opz, n.boot=10, size.boot=NULL, 
 %   jt=FALSE, nonParam=TRUE, random=FALSE, break.boot=n.boot)
%step.Ar.fit(obj, XREG, Z, PSI, opz, return.all.sol=FALSE)
%step.Ar.fit.boot(obj, XREG, Z, PSI, opz, n.boot=10, size.boot=NULL, jt=FALSE,
 %   nonParam=TRUE, random=FALSE, break.boot=n.boot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
vector of observations of length \code{n}.
}
  \item{x.lin, XREG}{
design matrix for standard linear terms.
}
  \item{Xtrue, Z}{
appropriate matrix including the stepmented variables whose breakpoints have to be estimated.
}
  \item{PSI}{
  appropriate matrix including the starting values of the breakpoints to be estimated.
}
  \item{ww,w}{
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
  \item{break.boot}{
  Integer, less than \code{n.boot}. If \code{break.boot} consecutive bootstrap samples lead to the same objective function, the algorithm stops without performing all \code{n.boot} 'trials'.
  This can save computational time considerably.
}
  \item{size.boot}{
  the size of the bootstrap resamples. If \code{NULL} (default), it is taken equal to the sample size.
  values smaller than the sample size are expected to increase perturbation in the bootstrap resamples.
}
  \item{jt}{
  logical. If \code{TRUE} the values of the stepmented variable(s) are jittered before fitting the model to the
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
%  \item{obj}{
%  the starting regression model where the stepmented relationships have to be added.
%}
%  \item{mfExt}{
%  the model frame.
%}

}
\details{
The functions call iteratively \code{lm.wfit} (or \code{glm.fit}) with proper design matrix depending on 
\code{XREG}, \code{Z} and \code{PSI}. \code{step.lm.fit.boot} (and \code{step.glm.fit.boot}) implements the bootstrap restarting idea discussed in
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
\code{\link{stepmented.lm}} or \code{\link{stepmented.glm}}
}
\examples{
##See ?stepmented
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{nonlinear }
