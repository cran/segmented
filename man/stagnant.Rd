\name{stagnant}
\alias{stagnant}
\docType{data}
\title{Stagnant band height data}
\description{
  The \code{stagnant} data frame has 28 rows and 2 columns.
}
\usage{data(stagnant)}
\format{
  A data frame with 28 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{log of flow rate in g/cm sec.}
    \item{\code{y}}{log of band height in cm}
  }
}
\details{
  Bacon and Watts report that such data were obtained by R.A. Cook during his investigation
  of the behaviour of stagnant surface layer height in a controlled flow of water.
}
\source{
  Bacon D.W., Watts D.G. (1971) Estimating the transistion between two intersecting straight lines.
  \emph{Biometrika} \bold{58}: 525 -- 534.
  
  Originally from the PhD thesis by R.A. Cook
}
%\references{
%  PhD thesis by R.A. Cook
%}
\examples{
data(stagnant)
## plot(stagnant)
}
\keyword{datasets}
