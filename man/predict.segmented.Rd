\name{predict.segmented}
\alias{predict.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for segmented model fits
}
\description{
Returns predictions and optionally associated quantities (standard errors or confidence intervals) from a fitted segmented model object.
}
\usage{
\method{predict}{segmented}(object, newdata, se.fit=FALSE, interval=c("none","confidence", "prediction"), 
            type = c("link", "response"), na.action=na.omit, level=0.95, .coef=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted segmented model coming from \code{segmented.(g)lm} or \code{segreg}.
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
Predictions on the link or response scale? Only if \code{object} is a segmented glm.
}
  \item{na.action}{
How to deal with missing data,  \emph{if} \code{newdata} include them.
}
  \item{level}{
The confidence level.
}
\item{.coef}{ The regression parameter estimates. If unspecified (i.e. \code{NULL}), it is computed internally by \code{coef()}.}
  \item{\dots}{
further arguments.
}
}
\details{
Basically \code{predict.segmented} builds the right design matrix accounting for breakpoint and passes it 
to \code{predict.lm} or \code{predict.glm} depending on the actual model fit \code{object}.
}
\value{
\code{predict.segmented} produces a vector of predictions with possibly associated standard errors or confidence intervals. 
See \code{predict.lm} or  \code{predict.glm}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
%\note{
%If \code{type="terms"}, \code{predict.segmented} returns predictions for each component of the segmented term. 
%Namely if `my.x' is the segmented variable, predictions for `my.x', `U1.my.x' and `psi1.my.x' are returned. These are
%meaningless individually, however their sum provides the predictions for the segmented term.
%}

\section{Warning }{
For segmented glm fits with offset, \code{predict.segmented} returns the fitted values \emph{including} the offset. 
}

\seealso{
\code{\link{segreg}}, \code{\link{segmented}}, \code{\link{plot.segmented}}, \code{\link{broken.line}}, \code{\link{predict.lm}}, \code{\link{predict.glm}}
}
\examples{
n=10
x=seq(-3,3,l=n)
set.seed(1515)
y <- (x<0)*x/2 + 1 + rnorm(x,sd=0.15)
segm <- segmented(lm(y ~ x), ~ x, psi=0.5)
predict(segm,se.fit = TRUE)$se.fit

#wrong (smaller) st.errors (assuming known the breakpoint)
olm<-lm(y~x+pmax(x-segm$psi[,2],0))
predict(olm,se.fit = TRUE)$se.fit
}
% \dontrun{..}
% KEYWORDS - R documentation directory.
\keyword{models} 
\keyword{regression}
