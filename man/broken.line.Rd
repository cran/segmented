\name{broken.line}
\alias{broken.line}
\title{ Fitted values for segmented relationships}
\description{
  Given a segmented model (typically returned by a \code{segmented} method), \code{broken.line}
  computes the fitted values (and relevant standard errors) for each `segmented' relationship.
}
\usage{
broken.line(ogg, term = NULL, link = TRUE, interc=TRUE, se.fit=TRUE)
}
\arguments{
  \item{ogg}{ A fitted object of class segmented (returned by any \code{segmented} method). }
  \item{term}{ Three options. A list (whose name should be one of the segmented covariates) 
    including values for which segmented predictions should be computed. A character meaning 
    the name of any segmented covariate in the model. \code{NULL} if the model includes a single
    segmented covariate. }
  \item{link}{ Should the predictions be computed on the scale of the link function? Default to \code{TRUE}. }
  \item{interc}{ Should the model intercept be added? (provided it exists).}
  \item{se.fit}{ If \code{TRUE} also standard errors for predictions are returned.}
}
\details{
  If \code{term=NULL} or \code{term} is a valid segmented covariate name, 
  predictions for each segmented variable are the relevant fitted values from the model. If \code{term} 
  is a (correctly named) list with numerical values, predictions corresponding to such specified values 
  are computed. If \code{link=FALSE} and \code{ogg} inherits from the class "glm", predictions and standard 
  errors are returned on the response scale. The standard errors come from the Delta method.
  Argument \code{link} is ignored whether \code{ogg} does not inherit from the class "glm".
    }
\value{
  A 2-component (if \code{se.fit=TRUE}) list representing predictions and standard errors for the segmented covariate values.
    }
%\references{ ~put references to the literature/web site here ~ }
\author{ Vito M. R. Muggeo }
%\note{ 
%This function will be probably removed in the next versions. See \code{predict.segmented} instead.
%}
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{segmented}}, \code{\link{predict.segmented}}, \code{\link{plot.segmented}}}
\examples{
set.seed(1234)
z<-runif(100)
y<-rpois(100,exp(2+1.8*pmax(z-.6,0)))
o<-glm(y~z,family=poisson)
o.seg<-segmented(o,seg.Z=~z,psi=.5)
\dontrun{plot(z,y)}
\dontrun{points(z,broken.line(o.seg,link=FALSE)$fit,col=2)} #ok, but use plot.segmented()!
    }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
