\name{selgmented}
\alias{selgmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selecting number of breakpoints in segmented regression
}
\description{
This function selects the number of breakpoints of the segmented relationship according to the BIC criterion or sequential hypothesis testing.
}
\usage{
selgmented(olm, seg.Z, alpha = 0.05, type = c("score", "davies", "bic"), 
  control = seg.control(), return.fit = TRUE, bonferroni = FALSE, Kmax = 2, msg=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{olm}{
A starting \code{lm} or \code{glm} object
}
  \item{seg.Z}{
A one-side formula for the segmented variable. Only one term can be included, and it can be omitted if \code{olm} includes just 
one covariate.}
  \item{alpha}{
The fixed type I error probability.
}
  \item{type}{
Which criterion should be used? Currently via sequential hypothesis testing \code{score} and \code{davies}, only \code{Kmax=2} is allowed. Using the BIC allows to compare models with more than 2 breakpoints.
}
  \item{control}{
See \code{\link{seg.control}}.
}
  \item{return.fit}{
If \code{TRUE}, the fitted model (with the selected number of breakpoints according to \code{type}) is returned.
}
  \item{bonferroni}{
If \code{TRUE}, the Bonferroni correction is employed, i.e. \code{alpha/Kmax} is always taken as threshold value to reject or not. If \code{FALSE}, \code{alpha} is used in the second level of hypothesis testing. 
}
  \item{Kmax}{
The maximum number of breakpoints being tested. If \code{type='bic'} any integer value can be specified, otherwise at most \code{Kmax=2} breakpoints can be tested via the Score or Davies statistics.
}
  \item{msg}{
If \code{FALSE} the final fit is returned silently with the selected number of breakpoints, otherwise the message including information about the selection procedure (i.e. the BIC values) is printed.
}
}
\details{
The function uses properly the functions \code{segmented}, \code{pscore.test} or \code{davies.test} to select the 'right' number
of breakpoints.
}
\value{The returned object depends on argument \code{return.fit}. If \code{FALSE}, the returned object is a list with some information on the compared models (i.e. the BIC values), otherwise a classical segmented object with the component \code{selection.psi} including the aforemntioned information. See \code{\link{segmented}} for details.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Muggeo V (2020) Selecting number of breakpoints in segmented regression: implementation in the R package segmented
https://www.researchgate.net/publication/343737604
}
\author{
Vito Muggeo
}
\note{
This is a sperimental function. Please use that with caution.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{segmented}}, \code{\link{pscore.test}}, \code{\link{davies.test}}
}
\examples{
## from ?segmented
\dontrun{
set.seed(12)
xx<-1:100
zz<-runif(100)
yy<-2+1.5*pmax(xx-35,0)-1.5*pmax(xx-70,0)+15*pmax(zz-.5,0)+rnorm(100,0,2)
dati<-data.frame(x=xx,y=yy,z=zz)
out.lm<-lm(y~x,data=dati)

os<-selgmented(out.lm) ## selects number of breakpoints via the Score test

os <-selgmented(out.lm, Kmax=3, type="bic") #BIC-based selection
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
