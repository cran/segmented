\name{seg}
\alias{seg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specifying a segmented term in the segreg formula
}
\description{
Function used to define the segmented term within the segreg formula. The function simply passes relevant information to proper fitter functions.
}
\usage{
seg(x, psi = NA, npsi = 1, est = NA, R = NA, fixed.psi = NULL, by = NULL, f.x = I)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The segmented (numeric) covariate
}
  \item{psi}{
Possible starting value(s) for the breakpoint(s). If omitted, \code{npsi} has to be specified.
}
  \item{npsi}{
The number of breakpoints to estimate. Default \code{npsi=1}.
}
  \item{est}{
Possible vector (of length equal to \code{npsi+1}) of 0/1 to indicate whether the slopes have to be estimated of fixed to zero.
If \code{NA}, default, means all ones, namely every slope has to be estimated. Consecutive zeroes are not allowed.
}
  \item{R}{
Matrix to constrain the slopes. If provided, it overwrites the matrix coming from the specification of \code{est}
}
  \item{fixed.psi}{
Possible \emph{fixed} breakpoint values. These will be accounted for \emph{in addition} to the estimated ones; \code{slope} and \code{plot.segmented} will account for them correctly.
}
  \item{by}{
A possible factor meaning an interaction with the segmented term. It means that a different segmented relationship is fitted within   each category.
}
  \item{f.x}{
an optional function meaning a function to apply to the covariate before fitting
}
}
\details{
The function is used within \code{\link{segreg}} to 'build' information about the segmented relationships to fit. 
}
\value{
The function simply returns the covariate with added attributes relevant to segmented term
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
\note{
If any value is provided in \code{fix.psi}, the corresponding slope difference coefficient will be labelled by \code{*.fixed.*}. The \code{\link{slope}} function will compute the 'right' slopes also accounting for the fixed breakpoints.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{segreg}}
}
\examples{
##see ?segreg
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
