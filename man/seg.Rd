\name{seg}
\alias{seg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specifying a segmented/stepmented term in the segreg/stepreg formula
}
\description{
Function used to define a segmented (stepmented) term within the segreg (stepreg) formula. The function simply passes relevant information to proper fitter functions.
}
\usage{
seg(x, npsi = 1, psi = NA, est = NA, R = NA, fixed.psi = NULL, by = NULL, f.x = I)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The segmented/stepmented (numeric) covariate
}
  \item{npsi}{
The number of breakpoints/jumpoints to estimate. Default to \code{npsi=1}. If \code{by} has been specified, the same \code{npsi} applies to all categories of the factor \code{by}; otherwise it can be vector, wherein the entries represent  the number of breakpoints of the segmented relationships within the categories of \code{by}. If \code{npsi} is specified as a vector, the corresponding \code{est} should be a list. The \code{npsi} starting values are computed according the specification of \code{quant} in \code{\link{seg.control}}.
}
  \item{psi}{
Numerical vector indicating possible starting value(s) for the breakpoint(s). When provided, \code{psi} overwrites \code{npsi}. If \code{by} has been specified, \code{psi} can be a list, wherein the components represent  the starting values of the segmented/stepmented relationships within the categories of \code{by}.
}
  \item{est}{
Possible vector (of length equal to \code{npsi+1}) of 0's and 1's to indicate whether the slopes have to be estimated or fixed to zero. \code{NA}, the default, means all ones, namely every slope is estimated. Consecutive zeroes are not allowed. If \code{by} has been specified, the same \code{est} applies to all categories of the factor \code{by}; otherwise \code{est} can be a list, wherein the components represent the slope constraints of the segmented relationships within the categories of \code{by}. Note that, when a grouping variable has been specified in \code{by}, \code{npsi} can be specified as a vector, and \code{est} should be a list.
}
  \item{R}{
Matrix to constrain the slopes. If provided, it overwrites the matrix (which is built internally) coming from the specification of \code{est}.
}
  \item{fixed.psi}{
Possible \emph{fixed} breakpoint values to be accounted for \emph{in addition} to the estimated ones; \code{slope} and \code{plot.segmented} will account for them correctly.
}
  \item{by}{
A possible \emph{factor} meaning an interaction with the segmented term \code{x}. Hence, if specified, a different segmented relationship is fitted within each category of \code{by}.
}
  \item{f.x}{
an optional function meaning a function to apply to the covariate before fitting
}
}
\details{
The function is used within \code{\link{segreg}} and \code{stepreg} to 'build' information about the segmented relationships to fit. 
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
