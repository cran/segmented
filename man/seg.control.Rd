\name{seg.control}
\alias{seg.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary for controlling segmented model fitting }
\description{
  Auxiliary function as user interface for 'segmented' fitting. Typically
     only used when calling any 'segmented' method (\code{segmented.lm} or \code{segmented.glm}).
}
\usage{
seg.control(toll = 1e-04, it.max = 20, display = FALSE, last = TRUE,
    maxit.glm = 25, h = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{toll}{ positive convergence tolerance. }
  \item{it.max}{ integer giving the maximal number of iterations. }
  \item{display}{ logical indicating if output should be produced for each iteration. }
  \item{last}{ logical indicating if output should include only the last fitted model.}
  \item{maxit.glm}{ integer giving the maximum number of inner IWLS iterations (see details). }
  \item{h}{ positive factor (from zero to one) modifying the increments in breakpoint estimation (see details). }
}
\details{
  Fitting a `segmented' model is attained via fitting iteratively standard GLMs. The number of iteration is 
  governed by \code{it.max}, while the (maximum) number of (inner) iterations to fit the GLM at each outer
  iteration is fixed via \code{maxit.glm}. Usually three-four inner iterations may be sufficient.\cr

  If \code{last=TRUE}, the object resulting from \code{segmented.lm} (or \code{segmented.glm}) is a
  list of fitted GLM; the i-th model is the segmented model with the values of the breakpoints at the i-th iteration. \cr

  Sometimes to stabilize the procedure, it can be useful to set \code{h<1} to reduce the increments in the breakpoint estimation. At each
  iteration the updated estimate is usually given by \code{psi.new=psi.old+increm}. By setting \code{h<1}
  (actually \code{min(abs(h),1)} is considered) causes the following alterations to the algorithm:
  (i) the actual maximum number of iterations is increased up to \code{it.max+round(it.max/2)} and
  (ii) the breakpoint update is computed via \code{psi.new=psi.old+h*increm}
  after the \code{it.max}th iteration.
}
\value{
  A list with the arguments as components.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

     #decrease the maximum number inner iterations and display the 
     #evolution of the (outer) iterations
     seg.control(display = TRUE, maxit.glm=4)
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
