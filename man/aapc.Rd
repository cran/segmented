\name{aapc}
\alias{aapc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average annual per cent change in segmented trend analysis
}
\description{
Computes the average annual per cent change  to summarize piecewise linear relationships in segmented regression models.
}
\usage{
aapc(ogg, parm, exp.it = FALSE, conf.level = 0.95, wrong.se = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ogg}{ the fitted model returned by \code{segmented}. }
  \item{parm}{ the \emph{single} segmented variable of interest. It can be missing if the model includes a single segmented covariate. If missing and      \code{ogg} includes several segmented variables, the first one is considered.}
  \item{exp.it}{logical. If \code{TRUE}, the per cent change is computed, namely \eqn{\exp(\hat\mu)-1}{exp(mu)-1} where 
  \eqn{\mu=\sum_j \beta_jw_j}{mu=\sum j bjwj}, see `Details'.}
  \item{conf.level}{the confidence level desidered.}
  \item{wrong.se}{logical, if \code{TRUE}, the `wrong'' standard error (as discussed in Clegg et al. (2009)) ignoring 
  uncertainty in the breakpoint estimate is returned as an attribute \code{"wrong.se"}.}
  \item{...}{further arguments to be passed on to \code{vcov.segmented()}, such as \code{var.diff} or \code{is}.}

}
\details{
 To summarize the fitted piecewise linear relationship, Clegg et al. (2009) proposed the 'average annual per cent change' (AAPC)
 computed as the sum of the slopes (\eqn{\beta_j}{beta_j}) weighted by corresponding covariate sub-interval width (\eqn{w_j}{w_j}), namely 
 \eqn{\mu=\sum_j \beta_jw_j}{mu=sum_j beta_j w_j}. Since the weights are the breakpoint differences, the standard error of the AAPC should account
 for uncertainty in the breakpoint estimate, as discussed in Muggeo (2010) and implemented by \code{aapc()}.
  }
\value{
  \code{aapc} returns a numeric vector including point estimate, standard error and confidence interval for the AAPC relevant to variable specified in \code{parm}.
}
\references{
     Clegg LX, Hankey BF, Tiwari R, Feuer EJ, Edwards BK (2009) Estimating average annual per cent change in trend analysis. 
     \emph{Statistics in Medicine}, \bold{28}; 3670-3682
    
    Muggeo, V.M.R. (2010) Comment on `Estimating average annual per cent change in trend analysis' by Clegg et al., 
    Statistics in Medicine; 28, 3670-3682. \emph{Statistics in Medicine}, \bold{29}, 1958--1960.
    }
\author{Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} }
\note{
\code{exp.it=TRUE} would be appropriate only if the response variable is the log of (any) counts.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
set.seed(12)
x<-1:20
y<-2-.5*x+.7*pmax(x-9,0)-.8*pmax(x-15,0)+rnorm(20)*.3
o<-lm(y~x)
os<-segmented(o, psi=c(5,12))
aapc(os)
}
\keyword{ regression }
