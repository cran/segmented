\name{slope.segmented}
\alias{slope.segmented}
\title{ Summary for slopes of segmented relationships }
\description{
  Computes summary of the slopes of each segmented relationship in the fitted model.
}
\usage{
slope.segmented(ogg, level = 0.95)
}
\arguments{
  \item{ogg}{ an object of class \code{"segmented"}, returned by any \code{segmented} method. }
  \item{level}{ the confidence level required. }
}
\details{
  To fit broken-line relationships, \code{segmented} uses a parameterization whose coefficients are not 
  the slopes. Therefore given an object \code{"segmented"}, \code{slope.segmented} computes point estimates,
  standard errors, t-values and confidence intervals of the slopes of each segmented 
  relationship in the fitted model.
}
\value{
  \code{slope.segmented} returns a list for each variable having a piecewise-linear relation 
  with the response in the fitted model. Each list is a matrix with number of rows equal to number of
  segments and five columns summarizing the results.
}
\references{ Muggeo, V.M.R. (2003) Estimating regression models with unknown break-points. 
    \emph{Statistics in Medicine} \bold{22}, 3055--3071.}

\author{ Vito M. R. Muggeo, \email{vito.muggeo@giustizia.it} }

\note{ The returned summary is based on limiting Gaussian distribution for the model parameters involved 
    in the computations. Sometimes even with large sample sizes such approximations are questionable 
    (e.g., with small difference-in-slopes parameters) and the results returned by \code{slope.segmented} 
    might be unreliable. Therefore is responsability of the user to gauge the applicability of such asymptotic 
    approximations. Anyway, the t values may be not assumed for testing purposes 
    and they should be used just as guidelines to assess the estimates' uncertainty.}

\seealso{ \code{\link{segmented}} }
\examples{
## A segmented relationship with three breakpoints
## 
x<-1:100
y<-2+1.5*pmax(x-35,0)-1.5*pmax(x-70,0)+rnorm(100,0,5)
dati<-data.frame(x=x,y=y)
out.lm<-lm(y~x,data=dati)
out.seg<-segmented(out.lm,Z=x,psi=c(20,80))
## the slopes of the three segments....
slope.segmented(out.seg)
rm(x,y,out.lm,out.seg)
}
\keyword{ regression }
