\name{down}
\alias{down}
\docType{data}
\title{ Down syndrome in babies}
\description{
The \code{down} data frame has 30 rows and 3 columns.
Variable \code{cases} means the number of babies with Down syndrome out of total number of births 
\code{births} for mothers with mean age \code{age}.
}
\usage{data(down)}
\format{
  A data frame with 30 observations on the following 3 variables.
  \describe{
    \item{\code{age}}{the mothers' mean age.}
    \item{\code{births}}{count of total births.}
    \item{\code{cases}}{count of babies with Down syndrome.}
  }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\source{
  Davison, A.C. and Hinkley, D. V. (1997) \emph{Bootstrap Methods
  and their Application}. Cambridge University Press.
}
\references{
  Geyer, C. J. (1991) Constrained maximum likelihood exemplified by isotonic 
  convex logistic regression. \emph{Journal of the American Statistical Association}
  \bold{86}, 717--724.
}
\examples{
data(down)
    }
\keyword{datasets}
