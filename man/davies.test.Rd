\name{davies.test}
\alias{davies.test}
\title{ Testing for a change in the slope }
\description{
  Given a generalized linear model, the Davies' test can be employed to test for a non-constant regression parameter
  in the linear predictor.
}
\usage{
davies.test(obj, seg.Z, k = 10, alternative = c("two.sided", "less", "greater"), 
    type=c("lrt","wald"), values=NULL, dispersion=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ a fitted model returned by \code{glm} or \code{lm}.  }
  \item{seg.Z}{ a formula with no response variable, such as \code{seg.Z=~x1}, indicating the 
      (continuous) segmented variable being tested. Only a single variable may be tested and a
      warning is printed when \code{seg.Z} includes two or more terms. }
  \item{k}{ number of points where the test should be evaluated. See Details. }
  \item{alternative}{ a character string specifying the alternative hypothesis. }
  \item{type}{ the test statistic to be used (only for GLM, default to lrt. 
  Ignored if \code{obj} is a simple linear model.}
  \item{values}{ optional. The evaluation points where the Davies approximation is computed. See Details for default values.}
  \item{dispersion}{ the dispersion parameter for the family to be used to compute the test statistic.
      When \code{NULL} (the default), it is inferred from \code{obj}. Namely it is taken as \code{1} for the
     Binomial and Poisson families, and otherwise estimated by the residual Chi-squared statistic (calculated from cases with
     non-zero weights) divided by the residual degrees of freedom.
       }
}
\details{
  \code{davies.test} tests for a non-zero difference-in-slope parameter of a segmented
  relationship. Namely, the null hypothesis is \eqn{H_0:\beta=0}{H_0:beta=0}, where \eqn{\beta}{beta} is the difference-in-slopes, 
  i.e. the coefficient of the segmented function \eqn{\beta(x-\psi)_+}{beta*(x-psi)_+}. The hypothesis of interest 
  \eqn{\beta=0}{beta=0} means no breakpoint. 
  Roughtly speaking, the procedure computes \code{k} `naive' (i.e. assuming
  fixed and known the breakpoint) test statistics for the difference-in-slope,
  seeks the `best' value and corresponding naive p-value (according to the alternative hypothesis), and then corrects 
  the selected (minimum) p-value by means of the \code{k} values of the test statistic. 
  If \code{obj} is a LM, the Davies (2002) test is implemented. This approach works even for small samples. 
  If \code{obj} represents a GLM fit, relevant methods are described in Davies (1987), and the Wald or the Likelihood ratio 
  test statistics can be used, see argument \code{type}. This is an asymptotic test.
#  The \code{k} evaluation points are \code{k} equally spaced values between the minimum and the maximum (excluded) values of the variable reported in \code{seg.Z}.
  The \code{k} evaluation points are \code{k} equally spaced values between the second and the second-last values of the variable reported in \code{seg.Z}.
}
\value{
  A list with class '\code{htest}' containing the following components:
  \item{method}{title (character)}
  \item{data.name}{the regression model and the segmented variable being tested}
  \item{statistic }{the point within the range of the covariate in \code{seg.Z} at which the maximum 
  (or the minimum if \code{alternative="less"}) occurs}
  \item{parameter }{number of evaluation points}
  \item{p.value }{the adjusted p-value}
  \item{process}{a two-column matrix including the evaluation points and corresponding values of the test statistic}
}
\references{
Davies, R.B. (1987) Hypothesis testing when a nuisance parameter is present only under the alternative.
    \emph{Biometrika} \bold{74}, 33--43. 
    
Davies, R.B. (2002) Hypothesis testing when a nuisance parameter is present only under the alternative: 
    linear model case. \emph{Biometrika} \bold{89}, 484--489. 
    }
\author{ Vito M.R. Muggeo }
\note{
  Strictly speaking,
  the Davies test is not confined to the segmented regression; the procedure can be applied when a nuisance parameter
  vanishes under the null hypothesis. The test is slightly conservative, as the computed p-value is actually an upper
  bound.
#  results should change slightly with respect to previous versions  where the evaluation points were computed 
#  as \code{seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)}.
  
  }

\section{Warning }{
    The Davies test is not aimed at obtaining the estimate of the breakpoint.
    The Davies test is based on \code{k} evaluation points, thus the value returned in the \code{statistic} component
    (and printed as "'best' at") is the best among the \code{k} points, and typically it will differ from the maximum likelihood estimate
    returned by \code{segmented}. Use \code{\link{segmented}} if you are interested in the point estimate.
  If \code{obj} is \code{"glm"} object with gaussian family, \code{davies.test()} uses an approximate
  test resulting in a typically smaller p-values in small samples.   Therefore to test for a breakpoint in linear 
  models with small samples, it is suggested to use \code{davies.test()} with objects of class "lm".
  However if the sample size is large (n>300), the exact Davies (2002) upper bound cannot be computed (as it relies on \code{gamma()} function
  ) and the \emph{approximate} upper bound of Davies (1987) is returned.
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
