\name{simulData}
\alias{simulData}
\docType{data}
\title{ Simulated data}
\description{
  Simulated dataset to illustrate the capability of segmented
}
\usage{data(simulData)}
\format{
  A data frame with 200 observations on the following 3 variables.
  \describe{
    \item{\code{y}}{the response variable}
    \item{\code{x}}{the explanatory variable}
    \item{\code{group}}{the grouping variable}
  }
}
\details{
  This dataset is employed to illustrate the capabilities of \code{segmented}
  in fitting segmented relationships with several breakpoints and constraints on the
  slope parameters.
}
\examples{
data(simulData)
## plot(simulData)
}
\keyword{datasets}
