\name{summary.segmented}
\alias{summary.segmented}
\alias{print.summary.segmented}
\title{ Summarizing  model fits for segmented regression }
\description{
  summary method for class \code{segmented}.
}
\usage{
summary.segmented(object, short = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class "segmented" }
  \item{short}{ logical indicating if the `short' summary should be printed }
  \item{\dots}{ further arguments }
}
\details{
  If short=TRUE only coefficients of the segmented relationships are printed.
}
\value{
A list (similar to one returned by \code{segmented.lm} or \code{segmented.glm}) with additional components:
  \item{psi }{estimated break-points and relevant (approximate) standard errors}
  \item{Ttable }{estimates and standard errors of the model parameters. This is similar
  to the matrix \code{coefficients} returned by \code{summary.lm} or \code{summary.glm},
  but without the rows corresponding to the breakpoints. Even the p-values relevant to the
  difference-in-slope parameters have been replaced by NA, since they are meaningless in
  this case, see \code{\link{davies.test}}}.
  \item{gap}{estimated coefficients, standard errors and t-values for the `gap' variables}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Vito M.R. Muggeo }
\seealso{ \code{\link{print.segmented}}, \code{\link{davies.test}} }
\examples{
  # continues example from segmented()
  # summary(segmented.model,short=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
