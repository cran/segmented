\name{print.segmented.lme}
\alias{print.segmented.lme}
\alias{fixef.segmented.lme}
\alias{logLik.segmented.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print method for the segmented.lme class }
\description{
  Printing and extracting the most important features of a segmented mixed model.
}
\usage{
\method{print}{segmented.lme}(x, digits = max(3, getOption("digits") - 3), ...)

\method{fixef}{segmented.lme}(object, ...)

\method{logLik}{segmented.lme}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{segmented.lme} }
  \item{digits}{ number of digits to be printed }
  \item{object}{object of class \code{segmented} }
  %\item{include.psi}{logical. If \code{TRUE}, the breakpoints are returned along with the regression coefficients}
  \item{\dots}{ arguments passed to other functions }
}
\author{ Vito M.R. Muggeo }

\seealso{ \code{\link{segmented.lme}}, \code{\link{summary.segmented.lme}} }

\keyword{ models }
