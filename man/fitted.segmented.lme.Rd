\name{fitted.segmented.lme}
\alias{fitted.segmented.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Fitted values for segmented mixed fits
}
\description{
Computes fitted values at different levels of nesting for segmented mixed objects
}
\usage{
\method{fitted}{segmented.lme}(object, level = 1, sort=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class \code{"segmented.lme"}
}
  \item{level}{
the level to be considered. Currently only levels 0 or 1 are allowed.
}
  \item{sort}{
If \code{TRUE}, the fitted values are sorted by the names of the 'id' levels.
}

  \item{\dots}{
Ignored
}
}
\details{
Currently it works only if \code{level=1}
}
\value{
A numeric object including the fitted values at the specified level of nesting.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
%\note{
%%  ~~further notes~~
%}

\section{Warning }{
All the functions for segmented mixed models (*.segmented.lme) are still at an experimental stage 
} 

\seealso{
\code{\link{summary.segmented.lme}}
}
%\examples{
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
