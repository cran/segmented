\name{glucoseData}
\alias{glucoseData}
\docType{data}
\title{
Response to Glucose Challenge
}
\description{
A glucose challenge is administered to individuals who may have diabetes to assess their ability to process glucose
}
\usage{data("glucoseData")}
\format{
  A data frame with 264 observations on the following 4 variables.
  \describe{
    \item{\code{subj}}{the subject id.}
    \item{\code{group}}{The group ('obese' vs 'control')}
    \item{\code{phos}}{The plasma inorganic phosphate values (mg/dl)}
    \item{\code{time}}{The observational times (hours after challenge).}
  }
}
\details{
Data that were obtained in an experiment with 20 obese and 13 control patients who took a standard glucose tolerance test. The response measurements are plasma inorganic phosphate (inversely related to glucose level) taken at different times.
}
\source{
Zerbe, G. O. (1979). Randomization analysis of the completely randomized design extended to growth and response curves. Journal of the American Statistical  Association, 74, 215-221.
}
\references{
Cudeck, R.,  Klebe, K. J. (2002). Multiphase mixed-effects models for repeated measures data. Psychological Methods, 7, 41-63.
}
\examples{
data(glucoseData)
## maybe str(glucoseData) ; plot(glucoseData) ...
}
\keyword{datasets}
