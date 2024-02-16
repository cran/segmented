\name{globTempAnom}
\alias{globTempAnom}
\docType{data}
\title{ Global temperature anomalies 1850-2023}
\description{
The \code{globTempAnom} data frame includes the global temperature anomalies from 1850 to 2023.
}
\usage{data(globTempAnom)}
\format{
  The included variables are (clearly).
  \describe{
    \item{\code{Year}}{the calendar year.}
    \item{\code{Anomaly}}{the temperature anomalies computed as differences with respect to the 20th century average (1901-2000) }
  }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\source{
  https://www.ncei.noaa.gov/access/monitoring/global-temperature-anomalies/anomalies
}
\references{
  There are several references using such dataset, e.g.

  Cahill, N., Rahmstorf, S., and Parnell, A. C. (2015). Change points of
global temperature. \emph{Environmental Research Letters}, 10: 1-6.
}
\examples{
data(globTempAnom)
    }
\keyword{datasets}
