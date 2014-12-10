\name{segmented-package}
\alias{segmented-package}
%\alias{segmented}
\docType{package}
\title{
Segmented relationships in regression models with breakpoints/changepoints estimation
}
\description{
Estimation of Regression Models with piecewise linear relationships having a
fixed number of break-points.
}
\details{
\tabular{ll}{
Package: \tab segmented\cr
Type: \tab Package\cr
Version: \tab 0.5-1.0\cr
Date: \tab 2014-12-10\cr
License: \tab GPL\cr
}
 Package \code{segmented} is aimed to estimate linear and generalized linear models (and virtually any regression model)
 having one or more segmented relationships in the linear predictor. Estimates of the slopes and
 of the possibly multiple breakpoints are provided. The package includes testing/estimating
 functions and methods to print, summarize and plot the results. \cr
 
 The algorithm used by \code{segmented} is \emph{not} grid-search. It is an iterative procedure (Muggeo, 2003) 
 that needs starting values \emph{only} for the breakpoint parameters and therefore it is quite efficient even
 with several breakpoints to be estimated. Moreover since version 0.2-9.0, \code{segmented} implements 
 the bootstrap restarting (Wood, 2001) to make the algorithm less sensitive to starting values. \cr

 Since version 0.5-0.0 a default method \code{segmented.dafault} has been added. It may be employed to include segmented relationships
 in \emph{general} regression models where specific methods do not exist. Examples include quantile and Cox regressions. See 
 examples in \code{\link{segmented.default}}.\cr
 
 A tentative approach to deal with unknown number of breakpoints
 is also provided, see option \code{stop.if.error} in \code{\link{seg.control}}.
}
\author{
Vito M.R. Muggeo <vito.muggeo@unipa.it>
  }
\references{

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

    }

\keyword{ regression }
\keyword{ nonlinear }
