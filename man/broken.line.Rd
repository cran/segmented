\name{broken.line}
\alias{broken.line}
\title{ Fitted values for segmented relationships}
\description{
  Given a segmented model (typically returned by a \code{segmented} method), \code{broken.line}
  computes the fitted values (and relevant standard errors) for the specified `segmented' relationship.
}
\usage{
broken.line(ogg, term = NULL, link = TRUE, interc=TRUE, se.fit=TRUE, isV=FALSE, 
  .vcov=NULL, .coef=NULL, ...)
}
\arguments{
  \item{ogg}{ A fitted object of class segmented (returned by any \code{segmented} method). }
  \item{term}{ Three options. i) A named list (whose name should be one of the segmented covariates in the model \code{ogg}) 
    including the covariate values for which segmented predictions should be computed; ii) a character meaning 
    the name of any segmented covariate in the model (and predictions corresponding to the observed covariate values are returned); 
    iii) It can be \code{NULL} if the model includes a single segmented covariate (and predictions corresponding to the observed covariate values are returned).}
  \item{link}{ Should the predictions be computed on the scale of the link function if \code{ogg} is a segmented glm fit? Default to \code{TRUE}. }
  \item{interc}{ Should the model intercept be added? (provided it exists).}
  \item{se.fit}{ If \code{TRUE} also standard errors for predictions are returned.}
  \item{isV}{ A couple of logicals indicating if the segmented terms \eqn{(x-\psi)_+}{(x-\psi)_+} and \eqn{I(x>\psi)}{I(x>\psi)} in the model matrix should be replaced by their smoothed counterparts when computing the standard errors. If a single logical is provided, it is applied to both terms.}
  \item{.vcov}{ Optional. The \emph{full} covariance matrix of estimates. If \code{NULL} (and \code{se.fit=TRUE}), the matrix is computed internally via \code{vcov.segmented()}.}
  \item{.coef}{ The regression parameter estimates. If unspecified (i.e. \code{NULL}), it is computed internally by \code{coef(ogg)}.}
  \item{...}{ Additional arguments to be passed on to \code{vcov.segmented()} when computing the standard errors for the predictions, namely
    \code{is}, \code{var.diff}, \code{p.df}. See \code{\link{summary.segmented}} and \code{\link{vcov.segmented}}.}
}
\details{
  If \code{term=NULL} or \code{term} is a valid segmented covariate name, 
  predictions for that segmented variable are the relevant fitted values from the model. If \code{term} 
  is a (correctly named) list with numerical values, predictions corresponding to such specified values 
  are computed. If \code{link=FALSE} and \code{ogg} inherits from the class "glm", predictions and possible standard 
  errors are returned on the response scale. The standard errors come from the Delta method.
  Argument \code{link} is ignored whether \code{ogg} does not inherit from the class "glm".
    }
\value{
  A list having one component if (if \code{se.fit=FALSE}), and two components (if \code{se.fit=TRUE}) list representing predictions and standard errors for the segmented covariate values.
    }
%\references{ ~put references to the literature/web site here ~ }
\author{ Vito M. R. Muggeo }
\note{ 
This function was written when there was not \code{predict.segmented} (which is more general).
}
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{segmented}}, \code{\link{predict.segmented}}, \code{\link{plot.segmented}}, \code{\link{vcov.segmented}}}
\examples{
set.seed(1234)
z<-runif(100)
y<-rpois(100,exp(2+1.8*pmax(z-.6,0)))
o<-glm(y~z,family=poisson)
o.seg<-segmented(o,seg.Z=~z)
\dontrun{plot(z,y)}
\dontrun{points(z,broken.line(o.seg,link=FALSE)$fit,col=2)} #ok, but use plot.segmented()!
    }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
