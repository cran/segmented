\name{selgmented}
\alias{selgmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Selecting the number of breakpoints in segmented regression
}
\description{
This function selects (and estimates) the number of breakpoints of the segmented relationship according to the BIC/AIC criterion or sequential hypothesis testing.
}
\usage{
selgmented(olm, seg.Z, Kmax=2, type = c("score", "bic", "davies", "aic"), 
  alpha = 0.05, control = seg.control(), refit = FALSE, stop.if = 5, 
  return.fit = TRUE, bonferroni = FALSE, msg = TRUE, plot.ic = FALSE, th = NULL, 
  G = 1, check.dslope = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{olm}{
A starting \code{lm} or \code{glm} object, or a simple numerical vector meaning the response variable.
}
  \item{seg.Z}{
A one-side formula for the segmented variable. Only one term can be included, and it can be omitted if \code{olm} is a (g)lm fit including just one numeric covariate. Also it might be omitted, and will be taken as 1,2..., if \code{olm} includes a single numeric variable. 
}
  \item{Kmax}{
The maximum number of breakpoints being tested. If \code{type='bic'} or \code{type='aic'}, any integer value can be specified; otherwise at most \code{Kmax=2} breakpoints can be tested via the Score or Davies statistics.
}
  \item{type}{
Which criterion should be used? Options \code{"score"} and \code{"davies"} allow to carry out sequential hypothesis testing with no more than 2 breakpoints (\code{Kmax=2}). Alternatively, the number of breakpoints can be selected via the BIC (or AIC) with virtually no upper bound for \code{Kmax}.
}
  \item{alpha}{
The fixed type I error probability when sequential hypothesis testing is carried out (i.e. \code{type='score'} or \code{'davies'}). It is also used when \code{type='bic'} (or \code{type='aic'}) and \code{check.dslope=TRUE} to remove the breakpoints based on the slope diffence t-value.
}
  \item{control}{
See \code{\link{seg.control}}.
}
  \item{refit}{
If \code{TRUE}, the final selected model is re-fitted using arguments in \code{control}, typically with bootstrap restarting. Set \code{refit=FALSE} to speed up computation (and possibly accepting near-optimal estimates). It is always \code{TRUE} if \code{type='score'} or \code{type='davies'}.
}
  \item{stop.if}{
An integer. If, when trying models with an increasing (when \code{G=1}) or decreasing (when \code{G>1}) number of breakpoints, \code{stop.if} consecutive fits exhibit higher AIC/BIC values, the search is interrupted. Set a large number, larger then \code{Kmax} say, if you want to assess the fits for all breakpoints \code{0, 1, 2, ..., Kmax}. Ignored if \code{type='score'} or \code{type='davies'}.
}
  \item{return.fit}{
If \code{TRUE}, the fitted model (with the number of breakpoints selected according to \code{type}) is returned.
}
  \item{bonferroni}{
  If \code{TRUE}, the Bonferroni correction is employed, i.e. \code{alpha/Kmax} (rather than \code{alpha}) is always taken as threshold value to reject or not. If \code{FALSE}, \code{alpha} is used in the second level of hypothesis testing. It is also effective when \code{type="bic"} (or \code{'aic'}) and \code{check.dslope=TRUE}, see Details. 
}
  \item{msg}{
If \code{FALSE} the final fit is returned silently with the selected number of breakpoints, otherwise the messages about the selection procedure (i.e. the BIC values), and possible warnings are also printed.
}
  \item{plot.ic}{
If \code{TRUE} the information criterion values with respect to the number of breakpoints are plotted. Ignored if \code{type='score'} or \code{type='davies'} or \code{G>1}. Note that if \code{check.dslope=TRUE}, the final number of breakpoints could differ from the one selected by the BIC/AIC leading to an inconsistent plot of the information criterion, see Note below.   
}
  \item{th}{
When a large number of breakpoints is being tested, it could happen that 2 estimated breakpoints are too close each other, and only one can be retained. Thus if the difference between two breakpoints is less or equal to \code{th}, one (the first) breakpoint is removed. Of course, \code{th} depends on the \code{x} scale: Integers, like 5 or 10, are appropriate if the covariate is the observation index. Default (\code{NULL}) means \code{th=diff(range(x))/100}. Set \code{th=0} if you are willing to consider even breakpoints very clode each other. Ignored if \code{type='score'} or \code{type='davies'}.
}
  \item{G}{
Number of sub-intervals to consider to search for the breakpoints when \code{type='bic'} or \code{'aic'}. See Details.
}
  \item{check.dslope}{
Logical. Effective only if \code{type='bic'} or \code{'aic'}. After the optimal number of breakpoints has been selected (via AIC/BIC), should the \eqn{t}{t} values of the slope differences be checked? If \code{TRUE}, the breakpoints corresponding to slope differences with a 'low' \eqn{t}{t} values will be removed. Note the model is re-fitted at each removal and a new check is performed. Simulation evidence suggests that such strategy leads to better results. See Details.
}
%  \item{a}{
%An additional tuning parameter for the BIC. \eqn{a=1}{a=1} provides the classical BIC.
%}
}
\details{
The function uses properly the functions \code{segmented}, \code{pscore.test} or \code{davies.test} to select the 'optimal' number of breakpoints \code{0,1,...,Kmax}. If \code{type='bic'} or \code{'aic'}, the procedure stops if the last \code{stop.if} fits have increasing values of the information criterion. Moreover, a breakpoint is removed if too close to other, actually if the difference between two consecutive estimates is less then \code{th}. Finally, if \code{check.dslope=TRUE}, breakpoints whose corresponding slope difference estimate is `small' (i.e. \eqn{p}-value larger then \code{alpha} or \code{alpha/Kmax}) are also removed. 

When \eqn{G>1}{G>1} the dataset is split into \eqn{G}{G} groups, and the search is carried out separately within each group. This approach is fruitful when there are many breakpoints not evenly spaced in the covariate range and/or concentrated in some sub-intervals. \code{G=3} or \code{4} is recommended based on simulation evidence. 

Note \code{Kmax} is always tacitely reduced in order to have at least 1 residual df in the model with \code{Kmax} changepoints. Namely, if \eqn{n=20}, the maximal segmented model has \code{2*(Kmax + 1)} parameters, and therefore the largest \code{Kmax} allowed is 8.

When \code{type='score'} or \code{'davies'}, the function also returns the 'overall p-value' coming from combing the single p-values using the Fisher method. The pooled p-value, however, does not affect the final result which depends on the single p-values only. 

}

\note{
If \code{check.dslope=TRUE}, there is no guarantee that the final model has the lowest AIC/BIC. Namely the model with the best A/BIC could have `non-significant' slope differences which will be removed (with the corresponding breakpoints) by the final model. Hence the possible plot (obtained via \code{plot.ic=TRUE}) could be misleading. See Example 1 below.   
}

\value{
The returned object depends on argument \code{return.fit}. If \code{FALSE}, the returned object is a list with some information on the compared models (i.e. the BIC values), otherwise a classical \code{'segmented'} object (see \code{\link{segmented}} for details) with the component \code{selection.psi} including the A/BIC values and\cr
 - if \code{refit=TRUE}, \code{psi.no.refit} that represents the breakpoint values before the last fit (with boot restarting)\cr
 - if \code{G>1}, \code{cutvalues} including the cutoffs values used to split the data.

%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Muggeo V (2020) Selecting number of breakpoints in segmented regression: implementation in the R package segmented
https://www.researchgate.net/publication/343737604
}
\author{
Vito M. R. Muggeo
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{segmented}}, \code{\link{pscore.test}}, \code{\link{davies.test}}
}
\examples{

set.seed(12)
xx<-1:100
zz<-runif(100)
yy<-2+1.5*pmax(xx-35,0)-1.5*pmax(xx-70,0)+15*pmax(zz-.5,0)+rnorm(100,0,2)
dati<-data.frame(x=xx,y=yy,z=zz)
out.lm<-lm(y~x,data=dati)

os <-selgmented(out.lm) #selection (Kmax=2) via the Score test (default)

os <-selgmented(out.lm, type="bic", Kmax=3) #BIC-based selection

\dontrun{
########################################
#Example 1: selecting a large number of breakpoints

b <- c(-1,rep(c(1.5,-1.5),l=15))
psi <- seq(.1,.9,l=15)
n <- 2000
x <- 1:n/n
X <- cbind(x, outer(x,psi,function(x,y)pmax(x-y,0)))
mu <- drop(tcrossprod(X,t(b)))
set.seed(113)
y<- mu + rnorm(n)*.022 
par(mfrow=c(1,2))

#select number of breakpoints via the BIC (and plot it)
o<-selgmented(y, Kmax=20, type='bic', plot.ic=TRUE, check.dslope = FALSE) 
plot(o, res=TRUE, col=2, lwd=3)
points(o)
# select via the BIC + check on the slope differences (default)
o1 <-selgmented(y, Kmax=20, type='bic', plot.ic=TRUE) #check.dslope = TRUE by default
#note the plot of BIC is misleading.. But the number of psi is correct 
plot(o1, add=TRUE, col=3)
points(o1, col=3, pch=3)

##################################################
#Example 2: a large number of breakpoints not evenly spaced.  

b <- c(-1,rep(c(2,-2),l=10))
psi <- seq(.5,.9,l=10)
n <- 2000
x <- 1:n/n
X <- cbind(x, outer(x,psi,function(x,y)pmax(x-y,0)))
mu <- drop(tcrossprod(X,t(b)))
y<- mu + rnorm(n)*.02 

#run selgmented with G>1. G=3 or 4 recommended. 
#note G=1 does not return the right number of breaks  
o1 <-selgmented(y, type="bic", Kmax=20, G=4)
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
