\name{model.matrix.segmented}
\alias{model.matrix.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design matrix for segmented fits
}
\description{
This function builds the model matrix for \code{segmented} fits.
}
\usage{
\method{model.matrix}{segmented}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A segmented fit
}
  \item{\dots}{
additional arguments
}
}
\details{
\code{model.matrix.segmented}
}
\value{
The design matrix for a segmented regression model with the specified formula and data
}
\author{
Vito Muggeo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{model.matrix}}
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
