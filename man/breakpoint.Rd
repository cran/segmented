\name{breakpoint}
\alias{breakpoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Breakpoint value by covariate pattern
}
\description{
Computes point estimate, its standard error and confidence interval for breakpoint parameter modelled as a function of one or more covariates in a segmented mixed model 
}
\usage{
breakpoint(ogg, by = NULL, conf.level = 0.95, vcov. = NULL, use.link=TRUE, ...) %level = 0,
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ogg}{
an object of class "segmented.lme" 
}
  \item{by}{
It is a named list indicating covariate names and corresponding values affecting the breakpoint location. 
  For instance, \code{by=list(group="2",z2=.2)}, provided that the model has been fitted by specifying \code{group} and \code{z2} in 
  \code{z.psi}
}
  \item{conf.level}{
the specified confidence level
}
  \item{vcov.}{
The \emph{full} covariance matrix of estimates. If unspecified (i.e. \code{NULL}), the covariance matrix is computed internally by \code{vcov(ogg)}
}
  \item{use.link}{
     If \code{TRUE}, the changepoint value is returned on the scale used in the model fit, otherwise values always on the linear scale are returned. 
}
  \item{\dots}{
Optional arguments to be passed to \code{\link{vcov.segmented.lme}}
}
}
\details{
If the logit link is used for the changepoints, the returned standard error is computed via the Delta method.
}
\value{
A 1-row matrix with 4 columns
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{segmented.lme}}
}
%\examples{
%}
