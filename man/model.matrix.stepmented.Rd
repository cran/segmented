\name{model.matrix.stepmented}
\alias{model.matrix.stepmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design matrix for stepmented fits
}
\description{
This function builds the model matrix for \code{stepmented} fits.
}
\usage{
\method{model.matrix}{stepmented}(object, k=NULL, apprx=c("cdf","abs"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A stepmented fit
}
  \item{k}{
The (negative) exponent of the sample size to approximate the absolute value; see \code{\link{vcov.stepmented}}
}
  \item{apprx}{
The approximation of the absolute value.
}
  \item{\dots}{
additional arguments
}
}
\details{
\code{model.matrix.stepmented}
}
\value{
The design matrix for a stepmented regression model with the specified formula and data
}
\author{
Vito Muggeo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{model.matrix}}, \code{\link{vcov.stepmented}} 
}
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
