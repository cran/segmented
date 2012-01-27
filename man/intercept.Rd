\name{intercept}
\alias{intercept}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Intercept estimates from segmented relationships
}
\description{
Computes the intercepts of each `segmented' relationship in the fitted model.
}
\usage{
intercept(ogg, parm, gap = TRUE, rev.sgn = FALSE, var.diff=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ogg}{
     an object of class "segmented", returned by any \code{segmented} method.
}
  \item{parm}{
the segmented variable whose intercepts have to be computed. If missing all the segmented variables in the model are considered. 
}
  \item{gap}{
  logical. should the intercepts account for the (possible) gaps? 
}
 \item{rev.sgn}{vector of logicals. The length should be equal to the length of \code{parm}, but it is recycled otherwise.
  when \code{TRUE} it is assumed that the current \code{parm} is `minus' the actual segmented variable,
    therefore the sign is reversed before printing. This is useful when a null-constraint has been set on the last slope.
 }
  \item{var.diff}{Currently ignored as only point estimates are computed.
  %logical. If \code{var.diff=TRUE} and there is a single segmented variable, the computed standard errors 
  %  are based on a sandwich-type formula of the covariance matrix. See Details in \code{\link{summary.segmented}}.
    } 
 }

\details{
  A broken-line relationship means that a regression equation exists in the intervals 
  `\eqn{min(x)}{min(x)} to \eqn{\psi_1}{psi1}', `\eqn{\psi_1}{psi1} to \eqn{\psi_2}{psi2}', and so on.
  \code{intercept} computes point estimates of the intercepts of the different regression equations 
  for each segmented relationship in the fitted model.
}
\value{
\code{intercept} returns a list of one-column matrices. Each matrix represents a segmented relationship.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Vito M. R. Muggeo, \email{vito.muggeo@unipa.it}} 
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
See also \code{\link{slope}} to compute the slopes of the different regression equations 
  for each segmented relationship in the fitted model.
}
\examples{
## see ?slope
\dontrun{
intercept(out.seg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

