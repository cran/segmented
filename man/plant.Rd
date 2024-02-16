\name{plant}
\alias{plant}
\docType{data}
\title{ Plan organ dataset}
\description{
  The \code{plant} data frame has 103 rows and 3 columns.
}
\usage{data(plant)}
\format{
  A data frame with 103 observations on the following 3 variables:
  \describe{
    \item{\code{y}}{measurements of the plant organ.}
    \item{\code{time}}{times where measurements took place.}
    \item{\code{group}}{three attributes of the plant organ, \code{RKV}, \code{RKW}, \code{RWC}.}
  }
}
\details{
  Three attributes of a plant organ measured over time where biological reasoning indicates
  likelihood of multiple breakpoints.  The data are scaled to the maximum value
  for each attribute and all attributes are measured at each time.
}
\source{
  The data have been kindly provided by Dr Zongjian Yang
  at School of Land, Crop and Food Sciences, The University of Queensland, Brisbane, Australia.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
\dontrun{
data(plant)

lattice::xyplot(y~time,groups=group,auto.key=list(space="right"), data=plant)

o<-segreg(y~ 0+group+seg(time, by=group, npsi=2), data=plant)
summary(o)

par(mfrow=c(1,2))
plot(y~time, data=plant)
plot(o, term=1:3, add=TRUE, leg=NA, psi.lines=TRUE) #add the lines to the current plot

plot(o, term=1:3, col=3:5, res.col=3:5, res=TRUE, leg="bottomright")

}
}
\keyword{datasets}
