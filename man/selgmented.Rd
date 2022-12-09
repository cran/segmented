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
selgmented(olm, seg.Z, alpha = 0.05, type = c("score", "davies", "bic", "aic"), 
  control = seg.control(), refit=TRUE, stop.if=6, return.fit = TRUE, 
  bonferroni = FALSE, Kmax = 2, msg = TRUE, plot.ic = FALSE, th = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{olm}{
A starting \code{lm} or \code{glm} object or a simple numerical vector meaning the response variable.
}
  \item{seg.Z}{
A one-side formula for the segmented variable. Only one term can be included, and it can be omitted if \code{olm} includes just one covariate.}
  \item{alpha}{
The fixed type I error probability.
}
  \item{type}{
Which criterion should be used? Options \code{score} and \code{davies} allow to carry out sequential hypothesis testing with no more than 2 breakpoints (\code{Kmax=2}). Alternatively, the number of breakpoints can be selected via the BIC (or AIC) with virtually no upper bound for \code{Kmax}.
}
%  \item{break.type}{
%The kind of breakpoint being tested. \code{1} is for piecewise-linear relationships, \code{2} means piecewise-constant, i.e. a step-function, %relationships. See \code{\link{pscore.test}}. 
%}
  \item{control}{
See \code{\link{seg.control}}.
}
  \item{refit}{If \code{TRUE}, the final selected model is re-fitted using arguments in \code{control}, typically with bootstrap restarting. Set \code{refit=FALSE} to speed up computation (and possibly accepting near-optimal estimates). Ignored if \code{type='score'} or \code{type='davies'}. 
}
  \item{stop.if}{ An integer. If \code{stop.if} fits provide higher AIC/BIC values the search is interrupted. Ignored if \code{type='score'} or \code{type='davies'}.
}
  \item{return.fit}{
If \code{TRUE}, the fitted model (with the number of breakpoints selected according to \code{type}) is returned.
}
  \item{bonferroni}{
If \code{TRUE}, the Bonferroni correction is employed, i.e. \code{alpha/Kmax} is always taken as threshold value to reject or not. If \code{FALSE}, \code{alpha} is used in the second level of hypothesis testing. 
}
  \item{Kmax}{
The maximum number of breakpoints being tested. If \code{type='bic'} or \code{type='aic'}, any integer value can be specified, otherwise at most \code{Kmax=2} breakpoints can be tested via the Score or Davies statistics.
}
  \item{msg}{
If \code{FALSE} the final fit is returned silently with the selected number of breakpoints, otherwise the message including information about the selection procedure (i.e. the BIC values) is printed.
}
  \item{plot.ic}{
If \code{TRUE} the information criterion values with respect to the number of breakpoints are plotted. Ignored if \code{type='score'} or \code{type='davies'}. 
}
  \item{th}{
When a large number of breakpoints is being tested, it could happen that 2 estimated breakpoints are too close each other, and only one can be retained. Thus if the difference between two breakpoints is less or equal to \code{th}, one (the first) is deleted. Ignored if \code{type='score'} or \code{type='davies'}. Of course, \code{th} depends on the \code{x} scale: Integers, like 5 or 10, are appropriate if the covariate is the observation index. Default (\code{NULL}) means \code{th=diff(range(x))/100}.
}
}
\details{
The function uses properly the functions \code{segmented}, \code{pscore.test} or \code{davies.test} to select the 'optimal' number of breakpoints \code{0,1,...,Kmax}. If \code{type='bic'} or \code{'aic'}, the procedure stops if the last \code{stop.if} fits have increasing values of the information criterion. 
}
\value{The returned object depends on argument \code{return.fit}. If \code{FALSE}, the returned object is a list with some information on the compared models (i.e. the BIC values), otherwise a classical segmented object with the component \code{selection.psi} including the aforementioned information. See \code{\link{segmented}} for details.
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

########################################
#Selecting a large number of breakpoints

b<-c(-1,rep(c(1.5,-1.5),l=15))
psi<-seq(.1,.9,l=15)
n=2000
x<-1:n/n
X<-cbind(x, outer(x,psi,function(x,y)pmax(x-y,0)))
mu<-drop(X%*%b)
y<-mu + rnorm(n)*.02 
par(mfrow=c(1,2))
o<-selgmented(y, Kmax=20, type='bic', plot.ic=TRUE, refit=FALSE)
plot(o, res=T, col=2, lwd=3)

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
