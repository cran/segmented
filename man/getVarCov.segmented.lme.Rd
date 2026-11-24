\name{getVarCov.segmented.lme}
\alias{getVarCov.segmented.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract variance-covariance matrix
}
\description{
Extract the variance-covariance matrix from a fitted segmented mixed model
}
\usage{
\method{getVarCov}{segmented.lme}(obj, ..., type = c("random.effects", "conditional", "marginal"))
}
\arguments{
  \item{obj}{
A fitted segmented mixed model 
}
  \item{\dots}{
Arguments to be passed to \code{getVarCov.lme}, such as \code{individuals}
}
  \item{type}{
The type of variance-covariance matrix; \code{marginal} is not allowed.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A variance-covariance matrix or a list of variance-covariance matrices.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
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
