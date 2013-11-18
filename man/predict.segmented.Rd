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
\method{predict}{segmented}(object, newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted segmented model coming from \code{segmented.lm} or \code{segmented.glm}.
}
  \item{newdata}{
different choices: i) an optional data frame in which to look for variables with which to predict; 
ii) a data frame including values only for any segmented covariate in the model; iii) a (unnamed) numeric vector meaning 
the values of the single segmented covariate in the model; iv) if omitted, the fitted values are returned.
}
  \item{\dots}{
further arguments passed to \code{predict.lm} or \code{predict.glm}. Usually these 
are \code{se.fit}, or \code{interval} or \code{type}.
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
\note{
If \code{type="terms"}, \code{predict.segmented} returns predictions for each component of the segmented term. 
Namely if `my.x' is the segmented variable, predictions for `my.x', `U1.my.x' and `psi1.my.x' are returned. These are
meaningless of course.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predict.lm}}, \code{\link{predict.glm}}, \code{\link{plot.segmented}}
}
\examples{
n=10
x=seq(-3,3,l=n)
set.seed(1515)
y <- (x<0)*x/2 + 1 + rnorm(x,sd=0.15)
segm <- segmented(lm(y ~ x), ~ x, psi=0.5)
predict(segm,se.fit = TRUE)$se.fit

#wrong (lower) st.errors (assuming known the breakpoint)
olm<-lm(y~x+pmax(x-segm$psi[,2],0))
predict(olm,se.fit = TRUE)$se.fit
}
% \dontrun{..}
% KEYWORDS - R documentation directory.
\keyword{models} 
\keyword{regression}
