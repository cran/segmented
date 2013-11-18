\name{davies.test}
\alias{davies.test}
\title{ Testing for a change in the slope }
\description{
  Given a generalized linear model, the Davies' test can be employed to test for a non-constant regression parameter
  in the linear predictor.
}
\usage{
davies.test(obj, seg.Z, k = 10, alternative = c("two.sided", 
    "less", "greater"), beta0=0, dispersion=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ a fitted model returned by \code{glm} or \code{lm}.  }
  \item{seg.Z}{ a formula with no response variable, such as \code{seg.Z=~x1}, indicating the 
      (continuous) segmented variable being tested. Only a single variable may be tested and a
      warning is printed when \code{seg.Z} includes two or more terms. }
  \item{k}{ number of points where the test should be evaluated. See Details. }
  \item{alternative}{ a character string specifying the alternative hypothesis. }
  \item{beta0}{ the null value of the difference-in-slope; default to zero meaning no breakpoint. 
  Values different from zero are allowed, but see the Warning section. }
  \item{dispersion}{ the dispersion parameter for the family to be used to compute the Wald statistic.
      When \code{NULL} (the default), it is inferred from \code{obj}. Namely it is taken as \code{1} for the
     binomial and Poisson families, and otherwise estimated by the residual Chi-squared statistic (calculated from cases with
     non-zero weights) divided by the residual degrees of freedom.
       }
}
\details{
  \code{davies.test} tests for a non-zero difference-in-slope parameter of a segmented
  relationship. Namely, the null hypothesis is \eqn{H_0:\beta=\beta_0}{H_0:beta=beta0}, where \eqn{\beta}{beta} is the difference-in-slopes, 
  i.e. the coefficient of the segmented function \eqn{\beta(x-\psi)_+}{beta*(x-psi)_+}, and \eqn{\beta_0}{beta0} is the `null' value 
  specified via the argument \code{beta0}. Typically the hypothesis of interest assumes \eqn{\beta_0=0}{beta0=0} namely no breakpoint. 
  Roughtly speaking, the procedure computes \code{k} `naive' (i.e. assuming
  fixed and known the breakpoint) Wald statistics for the difference-in-slope,
  seeks the `best' value (according to the alternative hypothesis), and then corrects the selected
  (minimum) p-value. The \code{k} evaluation points are \code{k} equally spaced values between the second and the second-last 
  values of the variable reported in \code{seg.Z}.
}
\value{
  A list with class '\code{htest}' containing the following components:
  \item{method}{title (character)}
  \item{data.name}{the regression model and the segmented variable being tested}
  \item{statistic }{the point within the range of the covariate in \code{seg.Z} at which the maximum 
  (or the minimum if \code{alternative="less"}) occurs}
  \item{parameter }{number of evaluation points}
  \item{p.value }{the adjusted p-value}
  \item{process}{a two-column matrix including the evaluation points and corresponding values of the statistic}
}
\references{Davies, R.B. (1987) Hypothesis testing when a nuisance parameter is present only under the alternative.
    \emph{Biometrika} \bold{74}, 33--43. }
\author{ Vito M.R. Muggeo }
\note{
  Strictly speaking,
  the Davies test is not confined to the segmented regression; the procedure can be applied when a nuisance parameter
  vanishes under the null hypothesis. The test is slightly conservative, as the computed p-value is actually an upper
  bound.
  }

\section{Warning }{
    The Davies test is not aimed at obtaining the estimate of the breakpoint.
    The Davies test is based on \code{k} evaluation points, thus the value returned in the \code{statistic} component
    (and printed as "`Best' at") is the best among the \code{k} points, and typically it will differ from the maximum likelihood estimate
    returned by \code{segmented}. Use \code{\link{segmented}} if you are interested in the point estimate.
    Technically the Davies approach does not fit in contexts when \eqn{\beta_0\neq 0}, namely to test for a known difference-in-slope value.
    Therefore results when \code{beta0} differs from zero should be taken with caution.
    }
%%\section{Warning }{Currently \code{davies.test} does not work if the fitted model \code{ogg}
%%  does not include the segmented variable \code{term} being tested.}


\examples{
\dontrun{set.seed(20)
z<-runif(100)
x<-rnorm(100,2)
y<-2+10*pmax(z-.5,0)+rnorm(100,0,2)
o<-lm(y~z+x)

davies.test(o,~z)
davies.test(o,~x)
  }
}
\keyword{ htest }
