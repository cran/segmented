\name{predict.stepmented}
\alias{predict.stepmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for stepmented model fits
}
\description{
Returns predictions and optionally associated quantities (standard errors or confidence intervals) from a fitted stepmented model object.
}
\usage{
\method{predict}{stepmented}(object, newdata, se.fit=FALSE, interval=c("none","confidence", "prediction"), 
            type = c("link", "response"), na.action=na.omit, level=0.95, .coef=NULL, 
            .vcov=NULL, apprx.fit=c("none","cdf"), apprx.se=c("cdf","none"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted stepmented model coming from \code{stepmented.lm} or \code{stepmented.glm}.
}
  \item{newdata}{
An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
  \item{se.fit}{
Logical. Should the standard errors be returned?
}
  \item{interval}{
Which interval? See \code{\link{predict.lm}}
}
  \item{type}{
Predictions on the link or response scale? Only if \code{object} is a stepmented glm.
}
  \item{na.action}{
How to deal with missing data,  \emph{if} \code{newdata} include them.
}
  \item{level}{
The confidence level.
}
  \item{.coef}{ The regression parameter estimates. If unspecified (i.e. \code{NULL}), it is computed internally by \code{coef()}.}
  \item{.vcov}{ The estimate covariance matrix. If unspecified (i.e. \code{NULL}), it is computed internally by \code{vcov.stepmented()}.}
  \item{apprx.fit}{ The approximation of the \eqn{(x>\hat\psi)} used to compute the predictions/fitted values of the piece-wise relationships. }
  \item{apprx.se}{ The same abovementioned approximation to compute the standard error.}
  \item{\dots}{
further arguments, for instance \code{k} to be passed to \code{\link{vcov.stepmented}}.
}
}
\details{
Basically \code{predict.stepmented} builds the right design matrix accounting for breakpoint and passes it 
to \code{predict.lm} or \code{predict.glm} depending on the actual model fit \code{object}.
}
\value{
\code{predict.stepmented} produces a vector of predictions with possibly associated standard errors or confidence intervals. 
See \code{predict.lm}, \code{predict.glm}, or \code{\link{predict.segmented}}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
%\note{
%If \code{type="terms"}, \code{predict.stepmented} returns predictions for each component of the stepmented term. 
%Namely if `my.x' is the stepmented variable, predictions for `my.x', `U1.my.x' and `psi1.my.x' are returned. These are
%meaningless individually, however their sum provides the predictions for the stepmented term.
%}

\section{Warning }{
For stepmented glm fits with offset obtained starting from the model \code{glm(.., offset=..)}, \code{predict.stepmented} returns the fitted values without considering the offset. 
}

\seealso{
\code{\link{stepmented}}, \code{\link{plot.stepmented}}, \code{\link{predict.lm}}, \code{\link{predict.glm}}
}
\examples{
n=10
x=seq(-3,3,l=n)
set.seed(1515)
y <- (x<0)*x/2 + 1 + rnorm(x,sd=0.15)
segm <- segmented(lm(y ~ x), ~ x, psi=0.5)
predict(segm,se.fit = TRUE)$se.fit
}
% \dontrun{..}
% KEYWORDS - R documentation directory.
\keyword{models} 
\keyword{regression}
