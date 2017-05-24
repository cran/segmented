\name{pscore.test}
\alias{pscore.test}
\title{ Testing for existence of one breakpoint}
\description{
  Given a (generalized) linear model, the (pseudo) Score statistic tests for the existence of one breakpoint.
}
\usage{
pscore.test(obj, seg.Z, k = 10, alternative = c("two.sided", "less", "greater"), 
    values=NULL, dispersion=NULL, df.t=NULL, more.break=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ a fitted model typically returned by \code{glm} or \code{lm}. Even an object returned 
  by \code{segmented} can be set. Offset and weights are allowed.}
  \item{seg.Z}{ a formula with no response variable, such as \code{seg.Z=~x1}, indicating the 
      (continuous) segmented variable being tested. Only a single variable may be tested and an
      error is printed when \code{seg.Z} includes two or more terms. \code{seg.Z} can be omitted if i)\code{obj} is a segmented fit with a single segmented covariate (and that variable is taken), or ii)if it is a "lm" or "glm" fit with a single covariate (and that variable is taken).}
  \item{k}{ optional. Number of points used to compute the pseudo Score statistic. See Details. }
  \item{alternative}{ a character string specifying the alternative hypothesis. }
  \item{values}{ optional. The evaluation points where the Score test is computed. See Details for default values.}
  \item{dispersion}{ optional. the dispersion parameter for the family to be used to compute the test statistic.
      When \code{NULL} (the default), it is inferred from \code{obj}. Namely it is taken as \code{1} for the
     Binomial and Poisson families, and otherwise estimated by the residual Chi-squared statistic in the model \code{obj} (calculated from cases with
     non-zero weights divided by the residual degrees of freedom).}
  \item{df.t}{ optional. The degress-of-freedom used to compute the p-value. When \code{NULL}, the df extracted from \code{obj} are used.}
  \item{more.break}{ optional logical. If \code{obj} is a segmented fit, \code{more.break=FALSE} tests for the actual breakpoint in \code{obj}, 
  while \code{more.break=TRUE} tests for an \emph{additional} breakpoint. Ignored when \code{obj} is not a segmented fit.}

}
\details{
  \code{pscore.test} tests for a non-zero difference-in-slope parameter of a segmented
  relationship. Namely, the null hypothesis is \eqn{H_0:\beta=0}{H_0:beta=0}, where \eqn{\beta}{beta} is the difference-in-slopes, 
  i.e. the coefficient of the segmented function \eqn{\beta(x-\psi)_+}{beta*(x-psi)_+}. The hypothesis of interest 
  \eqn{\beta=0}{beta=0} means no breakpoint. Simulation studies have shown that such Score test is more powerful than the Davies test (see reference) when the alternative hypothesis is `one changepoint'.
  
  The \code{dispersion} value, if unspecified,  is taken from \code{obj}. If \code{obj} represents the fit under the null hypothesis (no changepoint), the dispersion parameter estimate will be usually larger, leading to a (potentially severe) loss of power.  
  
  The \code{k} evaluation points are \code{k} equally spaced values in the range of the segmented covariate. \code{k} should not be small. 
  Specific values can be set via \code{values}. However I have found no important difference due to number and location of the evaluation points, thus  default is \code{k=10} equally-spaced points. 
}
\value{
  A list with class '\code{htest}' containing the following components:
  \item{method}{title (character)}
  \item{data.name}{the regression model and the segmented variable being tested}
  \item{statistic }{the point within the range of the covariate in \code{seg.Z} at which the maximum 
  (or the minimum if \code{alternative="less"}) occurs}
  \item{parameter }{number of evaluation points}
  \item{p.value }{the p-value}
  \item{process}{the alternative hypothesis set}
}
\references{
Muggeo, V.M.R. (2016) Testing with a nuisance parameter present only under the alternative:
    a score-based approach with application to segmented modelling. 
    \emph{J of Statistical Computation and Simulation}, \bold{86}, 3059--3067. 
    }
\author{ Vito M.R. Muggeo }

\seealso{See also \code{\link{davies.test}}. }

\examples{
\dontrun{
set.seed(20)
z<-runif(100)
x<-rnorm(100,2)
y<-2+10*pmax(z-.5,0)+rnorm(100,0,3)

o<-lm(y~z+x)

#testing for one changepoint
#use the simple null fit
pscore.test(o,~z) #compare with davies.test(o,~z)..

#use the segmented fit
os<-segmented(o, ~z)
pscore.test(os,~z) #smaller p-value, as it uses the dispersion under the alternative (from 'os') 

#test for the 2nd breakpoint in the variable z
pscore.test(os,~z, more.break=TRUE) 

  }
}
\keyword{ htest }
