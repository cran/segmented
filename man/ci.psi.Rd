\name{ci.psi}
\alias{ci.psi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Interval estimate for the breakpoints}
\description{
  \code{ci.psi} returns confidence intervals for the breakpoints from a fitted segmented model.
}
\usage{
ci.psi(ogg, conf.level = 0.95, Fieller = FALSE, digits=max(3, getOption("digits") - 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ogg}{ a fitted object of class "segmented". }
  \item{conf.level}{ confidence level of the interval. Default to 0.95. }
  \item{Fieller}{ if TRUE, the Fieller's method is applied. Currently unimplemented. }
  \item{digits}{ controls the number of digits to print when printing the output. }
}
\details{
  Currently \code{ci.psi} computes confidence limits using the approximate standard errors of the estimated
  breakpoints as reported in the `psi' component of the list returned by any \code{segmented} method. The resulting 
  confidence intervals are based on the asymptotic Normal distribution of the breakpoint
  estimator which is reliable just for clear-cut kink relationships. See Details in \code{\link{segmented}}.
}

\value{

    A list of matrices. Each matrix includes point estimate and confidence limits of the breakpoint(s) for each
    segmented variable in the model.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

 \seealso{ see also \code{\link{segmented}} }


\examples{

set.seed(10)
x<-1:100
z<-runif(100)
y<-2+1.5*pmax(x-35,0)-1.5*pmax(x-70,0)+10*pmax(z-.5,0)+rnorm(100,0,2)
out.lm<-lm(y~x)
o<-segmented(out.lm,seg.Z=~x+z,psi=list(x=c(30,60),z=.4)) 
ci.psi(o)

    }

\keyword{ regression }

\keyword{ nonlinear }
