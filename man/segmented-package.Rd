\name{segmented-package}
\alias{segmented-package}
%\alias{segmented}
\docType{package}
\title{
Segmented relationships in regression models with breakpoints / changepoints estimation
}
\description{
Estimation and Inference of Regression Models with piecewise linear relationships having a
fixed number of break-points. The estimation method is described in Muggeo (2003) <doi: 10.1002/sim.1545>.
}
\details{
\tabular{ll}{
Package: \tab segmented\cr
Type: \tab Package\cr
Version: \tab 1.4-0\cr
Date: \tab 2022-01-28\cr
License: \tab GPL\cr
}
 Package \code{segmented} is aimed to estimate linear and generalized linear models (and virtually any regression model)
 having one or more segmented relationships in the linear predictor. Estimates of the slopes and
 breakpoints are provided along with standard errors. The package includes testing/estimating
 functions and methods to print, summarize and plot the results. \cr
 
 The algorithm used by \code{segmented} is \emph{not} grid-search. It is an iterative procedure (Muggeo, 2003) 
 that needs starting values \emph{only} for the breakpoint parameters and therefore it is quite efficient even
 with several breakpoints to be estimated. Moreover since version 0.2-9.0, \code{segmented} implements 
 the bootstrap restarting (Wood, 2001) to make the algorithm less sensitive to starting values. \cr

 Since version 0.5-0.0 a default method \code{segmented.default} has been added. It may be employed to include segmented relationships
 in \emph{general} regression models where specific methods do not exist. Examples include quantile, Cox, and lme regressions (where the random effects do not refer to the breakpoints!). See examples in \code{\link{segmented.default}}.\cr
 
 Since version 1.0-0 the estimating algorithm has been slight modified and it appears to be much stabler (in examples with noisy segmented relationhips and flat log likelihoods)
 then previous versions.\cr

 Hypothesis testing (about the existence of the breakpoint) and confidence intervals are performed via appropriate methods and functions.
 
 A tentative approach to deal with unknown number of breakpoints is also provided, see option \code{fix.npsi} in \code{\link{seg.control}}. Also, as version 1.3-0, the \code{\link{selgmented}} function has been introduced to select the number of breakpoints via the BIC or sequential hypothesis testing.
}
\author{
Vito M.R. Muggeo <vito.muggeo@unipa.it>
  }
\references{

Muggeo, V.M.R. (2017) Interval estimation for the breakpoint in segmented regression: a smoothed score-based approach. 
\emph{Australian & New Zealand Journal of Statistics} \bold{59}, 311--322.

Muggeo, V.M.R. (2016) Testing with a nuisance parameter present only under the alternative: a score-based approach with application to segmented modelling. \emph{J of Statistical Computation and Simulation} \bold{86}, 3059--3067.

Davies, R.B. (1987) Hypothesis testing when a nuisance parameter is present only under the alternative.
    \emph{Biometrika} \bold{74}, 33--43.

Seber, G.A.F. and Wild, C.J. (1989) \emph{Nonlinear Regression}. Wiley, New York.

Bacon D.W., Watts D.G. (1971) Estimating the transistion between two intersecting straight lines.
  \emph{Biometrika} \bold{58}: 525 -- 534.

Muggeo, V.M.R. (2003) Estimating regression models with unknown break-points.
    \emph{Statistics in Medicine} \bold{22}, 3055--3071.

Muggeo, V.M.R. (2008) Segmented: an R package to fit regression 
  models with broken-line relationships. \emph{R News} \bold{8/1}, 20--25.
  
Muggeo, V.M.R., Adelfio, G. (2011) Efficient change point detection in genomic sequences of 
    continuous measurements. \emph{Bioinformatics} \bold{27}, 161--166.

Wood, S. N. (2001) Minimizing model fitting objectives that contain spurious local minima
    by bootstrap restarting. \emph{Biometrics} \bold{57}, 240--244. 

Muggeo, V.M.R. (2010) Comment on `Estimating average annual per cent change in trend analysis' by Clegg et al., Statistics in Medicine; 28, 3670-3682. 
\emph{Statistics in Medicine}, \bold{29}, 1958--1960.
    }

\keyword{ regression }
\keyword{ nonlinear }
