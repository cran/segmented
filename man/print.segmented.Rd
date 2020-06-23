\name{print.segmented}
\alias{print.segmented}
\alias{coef.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print method for the segmented class }
\description{
  Printing the most important features and coefficients (including the breakpoints) of a segmented model.
}
\usage{
\method{print}{segmented}(x, digits = max(3, getOption("digits") - 3), ...)

\method{coef}{segmented}(object, include.psi=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{segmented} }
  \item{digits}{ number of digits to be printed }
  \item{object}{object of class \code{segmented} }
  \item{include.psi}{logical. If \code{TRUE}, the breakpoints are returned along with the regression coefficients}
  \item{\dots}{ arguments passed to other functions }
}
\author{ Vito M.R. Muggeo }

\seealso{ \code{\link{summary.segmented}}, \code{\link{print.summary.segmented}} }

\keyword{ models }
