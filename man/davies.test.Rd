\name{davies.test}
\alias{davies.test}
\title{ Testing for a change in the slope }
\description{
  Given a generalized linear model, the Davies' test can be employed to test for a non-constant regression parameter
  in the linear predictor.
}
\usage{
davies.test(ogg, term, k = 10,
  alternative = c("two.sided", "less", "greater"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ogg}{ a fitted model returned by \code{glm} or \code{lm}.  }
  \item{term}{ a character string to mean the segmented variable being tested. }
  \item{k}{ number of points where the test should be evaluated. See details. }
  \item{alternative}{ a character string specifying the alternative hypothesis. }
}
\details{
  \code{davies.test} tests for a non zero difference-in-slope parameter of a segmented
  relationship. Roughtly speaking, the procedure computes \code{k} `naive' (i.e. assuming
  fixed and known the breakpoint) Wald statistics for the difference-in-slope,
  seeks the `best' value (according to the alternative hypothesis), and then corrects the selected
  (minimum) p-value. The k evaluation points are the quantiles of the variable \code{term}.
}
\value{
  A list with class 'htest' containing the following components:
  \item{method}{title (character)}
  \item{data.name}{the regression model and the segmented variable being tested}
  \item{statistic }{the point at which the maximum (or the minimum if \code{alternative="less"}) occurs}
  \item{parameter }{number of evaluation points}
  \item{p.value }{the adjusted p-value}
}
\references{Davies, R.B. (1987) Hypothesis testing when a nuisance parameter is present only under the alternative.
    \emph{Biometrika} \bold{74}, 33--43. }
\author{ Vito M.R. Muggeo }
\note{
  Strictly speaking,
  the Davies test is not confined to the segmented regression; the procedure can be applied when a nuicance parameter
  vanishes under the null hypothesis. The test is slightly conservative, as the computed p-value is actually an upper
  bound.
  }

\section{Warning }{Currently \code{davies.test} does not work if the fitted model \code{ogg}
  has been built without the argument \code{data}.}


\examples{
set.seed(20)
z<-runif(100)
x<-rnorm(100,2)
y<-2+10*pmax(z-.5,0)+rnorm(100,0,2)
o<-lm(y~z+x)

davies.test(o,"z")
davies.test(o,"x")
}
\keyword{ htest }
