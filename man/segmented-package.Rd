\name{segmented-package}
\alias{segmented-package}
%\alias{segmented}
\docType{package}
\title{
Segmented Relationships in Regression Models with Breakpoints / Changepoints Estimation (with Possibly Random Effects) 
}
\description{
Estimation and inference of regression models with piecewise linear relationships, also known as segmented regression models, with a number of break-points fixed or to be `selected'. Random effects changepoints are also allowed since version 1.6-0, and since version 2.0-0 it is also possible to fit regression models with piecewise constant (or `stepmented') relationships.
}
\details{
\tabular{ll}{
Package: \tab segmented\cr
Type: \tab Package\cr
Version: \tab 2.0-3\cr
Date: \tab 2024-02-16\cr
License: \tab GPL\cr
}
 Package \code{segmented} aims to estimate linear and generalized linear models (and virtually any regression model) having one or more segmented or stepmented relationships in the linear predictor. Estimates of the slopes and
 breakpoints are provided along with standard errors. The package includes testing/estimating
 functions and methods to print, summarize and plot the results. 
 
 The algorithms used by \code{segmented} are \emph{not} grid-search. They are iterative procedures (Muggeo, 2003; Fasola et al., 2018) that need starting values \emph{only} for the breakpoint parameters and therefore they are quite efficient even
 with several breakpoints to be estimated. Moreover since version 0.2-9.0, \code{segmented} implements 
 the bootstrap restarting (Wood, 2001) to make the algorithms less sensitive to the starting values (which can be also omitted by the user) . 

 Since version 0.5-0.0 a default method \code{segmented.default} has been added. It may be employed to include segmented relationships in \emph{general} regression models where specific methods do not exist. Examples include quantile, Cox, and lme regressions where the random effects do not refer to the breakpoints; see \code{segmented.lme} to include random changepoints.  \code{\link{segmented.default}} includes some examples.
 
 Since version 1.0-0 the estimating algorithm has been slight modified and it appears to be much stabler (in examples with noisy segmented relationhips and flat log likelihoods)
 then previous versions.

 Hypothesis testing (about the existence of the breakpoint) and confidence intervals are performed via appropriate methods and functions.
 
 A tentative approach to deal with unknown number of breakpoints is also provided, see option \code{fix.npsi} in \code{\link{seg.control}}. Also, as version 1.3-0, the \code{\link{selgmented}} function has been introduced to select the number of breakpoints via the BIC or sequential hypothesis testing.
 
Since version 1.6-0, estimation of segmented mixed models has been introduced, see \code{\link{segmented.lme}} and related function. 
Since version 2.0-0, it is possible to fit segmented relationships with constraints on the slopes, see \code{\link{segreg}}. 

Finally, since 2.0-0, it is possible to fit (G)LM wherein one or more covariates  have a stepmented (i.e. a step-function like) relationship, see \code{\link{stepmented}}.
}
\author{
Vito M.R. Muggeo <vito.muggeo@unipa.it>
  }
\references{
Muggeo V.M.R., Atkins D.C., Gallop R.J., Dimidjian S. (2014) 
Segmented mixed models with random changepoints: a maximum likelihood approach 
with application to treatment for depression study. \emph{Statistical Modelling}, \bold{14}, 293-313.

Muggeo, V.M.R. (2017) Interval estimation for the breakpoint in segmented regression: a smoothed score-based approach. 
\emph{Australian & New Zealand Journal of Statistics}, \bold{59}, 311--322.

Fasola S, Muggeo V.M.R., Kuchenhoff, H. (2018) A heuristic, iterative algorithm for change-point detection in abrupt change models, \emph{Computational Statistics}, \bold{2}, 997--1015.

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
