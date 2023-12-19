\name{pwr.seg}
\alias{pwr.seg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Analysis in segmented regression
}
\description{
Given the appropriate input values, the function computes the power (sample size) corresponding to the specifed sample size (power). If a segmented fit object is provided, the power is computed taking the parameter estimates as input values.
}
\usage{
pwr.seg(oseg, pow, n, z = "1:n/n", psi, d, s, n.range = c(10,300), 
    X = NULL, break.type=c("break","jump"), alpha = 0.01, round.n = TRUE, 
    alternative = c("two.sided", "greater", "less"), msg = TRUE, ci.pow=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oseg}{
The fitted segmented object. If provided, the power is computed at the model parameter estimates, and all the remaining arguments but \code{alternative} and \code{alpha} are ignored. 
}
  \item{pow}{
The desired power level. If provided \code{n} has to be missing 
}
  \item{n}{
The fixed sample size. If provided \code{pow} has to be missing
}
  \item{z}{
The covariate understood to have a segmented effect. Default is \code{"1:n/n"}, i.e. equispaced values in (0,1). More generally a string indicating the quantile function having \code{p} and possible other numerical values as arguments. For istance \code{"qunif(p,0,1)"},  \code{"qnorm(p,2,5)"}, or \code{"qexp(p)"}. \code{"qunif(p,1,n)"} can be also specified, but attention should be paid to guarantee \code{psi} within the covariate range. Finally, it could be also a numerical vector meaning the actual covariate, but \code{pow} has to be missing. Namely if the covariate is supplied (and n is known), only the relevant power can be estimated.
}
  \item{psi}{
The breakpoint value within the covariate range
}
  \item{d}{
The slope difference
}
  \item{s}{
The response standard deviation
}
  \item{n.range}{
When \code{pow} is provided and the relevant sample size estimate has to be returned, the function evaluates 50 sample sizes equally spaced in \code{n.range}. However the function can also compute, via spline interpolation, sample sizes outside the   
specified range.
}
  \item{X}{
The design matrix including additional linear variables in the regression equation. Default to \code{NULL} which means intercept and linear term for the segmented covariate.
}
  \item{break.type}{
Type of breakpoint. \code{break.type='break'} means piecewise linear (segmented), \code{break.type='jump'} refers to piecewise constant. 
}
  \item{alpha}{
The type-I error probability. Default to 0.01.
}
  \item{round.n}{
logical. If \code{TRUE} the (possible) returned sample size value is rounded.
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" or "less". Note, this refers to the sign of the slope difference.
}
  \item{msg}{
logical. If \code{TRUE} the output is returned along with a simple message, otherwise only the values are returned 
}
  \item{ci.pow}{
Numerical. If \code{oseg} has been supplied, \code{ci.pow} replicates are drawn to build a 95\% confidence interval for the power. 
}
}
\details{
The function exploits the sampling distribution of the pseudo Score statistic under the alternative hypothesis of one breakpoint.
}
\value{
The computed power \emph{or} sample size, with or without message (depending on \code{msg})
}
\references{
D'Angelo N, Muggeo V.M.R. (2021) Power analysis in segmented regression, working paper \cr
  https://www.researchgate.net/publication/355885747.

Muggeo, V.M.R. (2016) Testing with a nuisance parameter present only under the alternative:
    a score-based approach with application to segmented modelling. 
    \emph{J of Statistical Computation and Simulation}, \bold{86}, 3059--3067.

}
\author{
Nicoletta D'Angelo and Vito Muggeo
}
\note{
Currently the function assumes just 1 breakpoint in one covariate
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pscore.test}}
}
\examples{
  ## pwr.seg(pow=.7, psi=.5, d=1.5, s=.5) #returns the sample size
  ## pwr.seg(n=219, psi=.5, d=1.5, s=.5) #returns the power
  ## pwr.seg(n=20,z="qnorm(p, 2,5)", psi=3, d=.5, s=2) #the covariate is N(2,5)
  ## pwr.seg(n=20,z="qexp(p)", psi=.1, d=.5, s=.1) #the covariate is Exp(1)

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
