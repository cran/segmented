\name{stepmented}
\alias{stepmented}
\alias{stepmented.lm}
\alias{stepmented.glm}
\alias{stepmented.ts}
\alias{stepmented.numeric}
%\alias{stepmented.default}
%\alias{stepmented.Arima}
%\alias{print.stepmented}
%\alias{summary.stepmented}
%\alias{print.summary.stepmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
stepmented relationships in regression models
}
\description{
  Fits regression models with stepmented (i.e. piecewise-constant) relationships between the response and one or more explanatory variables. Break-point estimates are provided.
}
\usage{
stepmented(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(), 
   keep.class=FALSE, var.psi=TRUE, ...)
%\method{stepmented}{default}(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(), model = TRUE, keep.class=FALSE, ...)
\method{stepmented}{lm}(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(),
   keep.class=FALSE, var.psi=TRUE, ...)

\method{stepmented}{glm}(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(),
    keep.class=FALSE, var.psi=TRUE, ...)

\method{stepmented}{numeric}(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(),
    keep.class=FALSE, var.psi=TRUE, ..., 
    pertV=0, centerX=FALSE, adjX=NULL, weights=NULL)

\method{stepmented}{ts}(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(),
    keep.class=FALSE, var.psi=TRUE, ..., 
    pertV=0, centerX=FALSE, adjX=NULL)


%\method{stepmented}{Arima}(obj, seg.Z, psi, npsi, fixed.psi=NULL, control = seg.control(),
%    model = TRUE, keep.class=FALSE, ...)
}
%- ------------------------------>>>>>  "Arima"
\arguments{
\item{obj}{A standard `linear' regression model of class "lm" or "glm". Alternatively, a simple "ts" object or a simple data vector may be supplied. 
}
\item{seg.Z}{ the stepmented variables(s), i.e. the numeric covariate(s) understood to have a piecewise-constant relationship with response. It is a formula with no response variable, such as \code{seg.Z=~x} or \code{seg.Z=~x1+x2}. Currently, formulas involving functions, 
such as \code{seg.Z=~log(x1)}, or selection operators, such as \code{seg.Z=~d[,"x1"]} or \code{seg.Z=~d$x1}, are \emph{not} allowed. Also, variable names formed by \verb{U} or \verb{V} only (with or without numbers ) are not permitted. If missing, the index variable \code{id=1,2,..,n} is used. For \code{stepmented.ts}, \code{seg.Z} is usually unspecified as the (time) covariate is obtained by the \code{ts} object itself.
}
\item{psi}{ starting values for the breakpoints to be estimated. If there is a single stepmented variable specified in \code{seg.Z}, \code{psi} can be a numeric vector, and it can be missing  when 1 breakpoint has to be estimated (and the median of the stepmented variable is used as a starting value). If \code{seg.Z} includes several covariates, \code{psi} has to be specified as a \emph{named} list of vectors whose names have to match the variables in the \code{seg.Z} argument. Each vector of such list includes starting values for the break-point(s) for the corresponding variable in \code{seg.Z}. A \code{NA} value means that `\code{K}' quantiles (or equally spaced values) are used as starting values; \code{K} is fixed via the \code{\link{seg.control}} auxiliary function. 
      }
\item{npsi}{A named vector or list meaning the \emph{number} (and not locations) of breakpoints to be estimated. The starting values will be internally computed via the quantiles or equally spaced values, as specified in argument \code{quant} in \code{\link{seg.control}}. \code{npsi} can be missing and \code{npsi=1} is assumed for all variables specified in \code{seg.Z}. If \code{psi} is provided, \code{npsi} is ignored.}
\item{fixed.psi}{An optional named list including the breakpoint values to be kept fixed during the estimation procedure. The names should be a subset of (or even the same) variables specified in \code{seg.Z}. If there is a single variable in \code{seg.Z}, a simple numeric vector can be specified. Note that, in addition to the values specified here, \code{stepmented} will estimate additional breakpoints. To keep fixed all breakpoints (to be specified in \code{psi}) use \code{it.max=0} in \code{\link{seg.control}}
}
\item{control}{ a list of parameters for controlling the fitting process.
      See the documentation for \code{\link{seg.control}} for details. 
      }
  %\item{model}{logical value indicating if the model.frame should be returned.}
  \item{keep.class}{logical value indicating if the final fit returned by \code{stepmented.default} should keep the class '\code{stepmented}' (along with the class of the original fit \code{obj}). Ignored by the stepmented methods. 
  }
  \item{\dots}{ optional arguments (to be ignored safely). Notice specific arguments relevant to the original call (via \code{lm} or \code{glm} for instance), such as \code{weights} or \code{offet}, have to be included in the starting model \code{obj}.
}
  \item{pertV}{
  Only for \code{stepmented.ts} and \code{stepmented.numeric}.
  }
  \item{centerX}{
  Only for \code{stepmented.ts} and \code{stepmented.numeric}. If \code{TRUE}, the covariate is centered before fitting.
  }
  \item{adjX}{
  Only for \code{stepmented.ts} and \code{stepmented.numeric}. If the response vector leads to covariate with large values (such as years for ts objects), \code{adjX=TRUE} will shift the covariate to have a zero origin. Default is \code{NULL} which means \code{TRUE} if the minimum of covariate is 1000 or larger.
}
\item{var.psi}{
logical. If \code{TRUE} (default), the estimate covariance matrix is also computed via \code{\link{vcov.stepmented}}, thus the breakpoint standard errors are also included in the \code{psi} component of the returned object.
}
  \item{weights}{
  possible weights to include in the estimation process (only for \code{stepmented.numeric}).
  }
}
\details{
  Given a linear regression model (usually of class "lm" or "glm"), stepmented tries to estimate
  a new regression model having piecewise-constant (i.e. step-function like) relationships with the variables specified in \code{seg.Z}.
  A \emph{stepmented} relationship is defined by the mean level
  parameters and the break-points where the mean level changes. The number of breakpoints
  of each stepmented relationship depends on the \code{psi} argument, where initial
  values for the break-points must be specified. The model
  is estimated simultaneously yielding point estimates and relevant approximate
   standard errors of all the model parameters, including the break-points.

  \code{stepmented} implements the algorithm described in Fasola et al. (2018) along with  bootstrap restarting (Wood, 2001) to escape local optima. The procedure turns out to be particularly appealing and probably efficient when there are two or more covariates exhibiting different change points to be estimated.

}
\value{
  The returned object depends on the \code{last} component returned by \code{seg.control}.
  If last=TRUE, the default, stepmented returns an object of class "stepmented" which inherits
  from the class "lm" or "glm" depending on the class of \code{obj}. Otherwise a list is
  returned, where the last component is the fitted model at the final iteration,
  see \code{\link{seg.control}}. \cr

An object of class "stepmented" is a list containing the components of the
original object \code{obj} with additionally the followings:
  \item{psi}{estimated break-points and relevant (approximate) standard errors}
  \item{it}{number of iterations employed}
  \item{epsilon}{difference in the objective function when the algorithm stops}
  \item{model}{the model frame}
  \item{psi.history}{a list or a vector including the breakpoint estimates at each step}
  \item{seed}{the integer vector containing the seed just before the bootstrap resampling. 
     Returned only if bootstrap restart is employed}
  \item{..}{Other components are not of direct interest of the user}
  }

%\section{ Warning }{
%It is well-known that the log-likelihood function for the 
%break-point may be not concave, especially 
%for poor clear-cut kink-relationships. In these circumstances, the initial guess
% for the break-point, i.e. the \code{psi} argument, must be provided with care. 
% For instance visual 
%inspection of a, possibly smoothed, scatter-plot is usually a good way to obtain some idea on breakpoint location. 
%However bootstrap restarting, implemented since version 0.2-9.0, is relatively more robust to starting values specified 
%in  \code{psi}. Alternatively an automatic procedure may be implemented by specifying \code{psi=NA} and \code{fix.npsi=FALSE} in %\code{\link{seg.control}}: experience suggests to increase the number of iterations via \code{it.max} in \code{seg.control()}. This automatic procedure, however, is expected to overestimate the number of breakpoints. 
%}

\note{
The component \code{psi.rounded} of the fit object includes the rounded changepoint values which are usually taken as the final estimates. \code{print.stepmented} and \code{print.sumamry.stepmented}and \code{confint.stepmented} return the rounded values.

Also:
\enumerate{

\item The algorithm will start if the \code{it.max} argument returned by \code{seg.control}
  is greater than zero. If \code{it.max=0} \code{stepmented} will estimate a new linear model with
 break-point(s) fixed at the values reported in \code{psi}.
 
\item In the returned fit object, `U.' is put before the name of the stepmented 
variable to indicate the difference in the mean levels. \code{\link{slope}} can be used to compute the actual mean levels corresponding to the different intervals.

\item Currently methods specific to the class \code{"stepmented"} are
    \itemize{
  \item \code{\link{print.stepmented}}
  \item \code{\link{summary.stepmented}}
  \item \code{\link{print.summary.stepmented}}
  \item \code{\link{plot.stepmented}}
  \item \code{\link{confint.stepmented}}
  \item \code{\link{vcov.stepmented}}
  \item \code{\link{lines.stepmented}}
  %\item \code{\link{predict.stepmented}}
  %\item \code{\link{points.stepmented}}
  %\item \code{\link{coef.stepmented}}
            }
Others are inherited from the class \code{"lm"} or \code{"glm"} depending on the
 class of \code{obj}.

     }
}


\references{
Fasola S, Muggeo VMR, Kuchenhoff H (2018) A heuristic, iterative algorithm for change-point detection in abrupt change models, \emph{Computational Statistics} \bold{33},  997--1015
  }
 
\author{ Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} (based on original code by Salvatore Fasola)}

\seealso{ \code{\link{segmented}} for segmented regression, \code{\link{lm}}, \code{\link{glm}} }

\examples{

n=20
x<-1:n/n
mu<- 2+ 1*(x>.6)
y<- mu + rnorm(n)*.8

#fitting via regression model
os <-stepmented(lm(y~1),~x)

y<-ts(y)
os1<- stepmented(y)  #the 'ts' method
os2<- stepmented(y, npsi=2)
#plot(y)
#plot(os1, add=TRUE)
#plot(os2, add=TRUE, col=3:5)


### Example with (poisson) GLM
y<- rpois(n,exp(mu))
o<-stepmented(glm(y~1,family=poisson))
plot(o, res=TRUE)

\dontrun{

### Example with (binary) GLM (example from the package stepR)
set.seed(1234)
y <- rbinom(200, 1, rep(c(0.1, 0.7, 0.3, 0.9), each=50))
o<-stepmented(glm(y~1,family=binomial), npsi=3)
plot(o, res=TRUE)

### Two stepmented covariates (with 1 and 2 breakpoints); z has also an additional linear effect
n=100
x<-1:n/n
z<-runif(n,2,5)
mu<- 2+ 1*(x>.6)-2*(z>3)+3*(z>4)+z
y<- mu + rnorm(n)*.8

os <-stepmented(lm(y~z),~x+z, npsi=c(x=1,z=2))
os
summary(os)

## see ?plot.stepmented
}
}

%# An example using the Arima method:
%\dontrun{
%n<-50
%idt <-1:n #the time index
%
%mu<-50-idt +1.5*pmax(idt-30,0)
%set.seed(6969)
%y<-mu+arima.sim(list(ar=.5),n)*3.5
%
%o<-arima(y, c(1,0,0), xreg=idt)
%os1<-stepmented(o, ~idt, control=seg.control(display=TRUE))
%
%#note using the .coef argument is mandatory!
%slope(os1, .coef=os1$coef)
%plot(y)
%plot(os1, add=TRUE, .coef=os1$coef, col=2)
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{nonlinear }
