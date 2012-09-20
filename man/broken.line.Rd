\name{broken.line}
\alias{broken.line}
\title{ Fitted values for segmented relationships}
\description{
  Given a segmented model (typically returned by a \code{segmented} method), \code{broken.line}
  computes the fitted values for each `segmented' relationship.
}
\usage{
broken.line(ogg, term = NULL, gap = FALSE, link = TRUE, 
    interc=TRUE)
}
\arguments{
  \item{ogg}{ A fitted object of class segmented (returned by any \code{segmented} method). }
  \item{term}{ A character meaning for which segmented variable prediction should be computed.}
  \item{gap}{ Should the `gaps' of the segmented relationships be plotted? Default to \code{FALSE}.}
  \item{link}{ Should the predictions be computed on the scale of the link function? Default to \code{TRUE}. }
  \item{interc}{ Should the model intercept be added? (provided it exists).}
}
\details{
  If \code{term=NULL} predictions for each segmented variable in the fitted model are computed. Argument \code{link}
  is ignored whether \code{ogg} does not inherit from the class "glm".
    }
\value{
  A matrix whose columns represent predictions for the segmented variables.
    }
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{segmented}}, \code{\link{predict.glm}}  }
\examples{
set.seed(1234)
z<-runif(100)
y<-rpois(100,exp(2+1.8*pmax(z-.6,0)))
o<-glm(y~z,family=poisson)
o.seg<-segmented(o,seg.Z=~z,psi=list(z=.5))
\dontrun{plot(z,y)}
\dontrun{points(z,broken.line(o.seg,link=FALSE),col=2,pch=20)}
    }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
