\name{sim.seg}
\alias{sim.seg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate data for (generalized) segmented (mixed) models
}
\description{
Function to simulate segmented mixed data, namely longitudinal data with individual profiles having  
piecewise linear relationships with a breakpoint. Heterogeneities (via multinormal random effects) 
are allowed for each model parameter: intercept, left and right slopes and breakpoint.
}
\usage{
sim.seg(N, ni, fixEff = c(2, -0.2, 0.3, 0.4), S = diag(c(0.05, 0.05, 0.05, 0.07)^2), 
  seg.Z=NULL, psii, fam=c("gauss","pois","bin"), z.psi, x.diff, 
  sd.eps = 0.01, seed, display = FALSE, eta0=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
number of clusters (usually subjects).
}
  \item{ni}{
number of measurements per subject. Three options. If \code{length(ni)=N}, it represents the actual number of measurements per each subject; 
if scalar (and larger than 1), exactly \code{ni} measurements per every subject are generated; otherwise the number of measurements 
for the \code{N} `subjects' is computed as \code{sample(min(n):max(n), size=N, replace=TRUE)}.
If \code{ni=1}, simple (i.e. without repeated measurements) segmented data are generated.
}
  \item{fixEff}{
A vector of four fixed effects: intercept, left slope, differences-in-slopes and changepoint. If the subject-specific changepoints are provided (via \code{psii}), \code{fixEff} should have 3 entries.
}
  \item{S}{
The covariance matrix of the random effects. The dimension should match the length of \code{fixEff}. Ignored if \code{ni=1}.
}
%  \item{psi.link}{
%A character to mean which link for the changepoint submodel should be used. Allowed values are 
%\code{"identity"} or \code{"logit"}.
%}
  \item{seg.Z}{
The numeric segmented covariate. Three options: i) a numeric vector; ii) a character to be evaluated to generate the variable; iii) \code{NULL} to generate uniform values in (0,1), differently among the subjects. If necessary, such vector is truncated or recycled to match the correct length. 
}

  \item{psii}{
An optional numeric vector (of length \code{N}) of subject-specific changepoints.
}
  \item{z.psi}{
Optional, to specify if the changepoint depends on some covariates. It is a two-components list including 
the coefficients and the covariate(s) affecting the changepoints. 
If there are p covariates, these should be organized in a \eqn{N \times p} matrix, otherwise a simple \eqn{N}-dimensional vector. 
\code{z.psi} can also be a single scalar representing the coefficient of the covariate obtained as \code{rep(c(0, 1), length = N)}. 
See the example below.
}
  \item{x.diff}{
Optional. Similar to \code{z.psi}, but for the slope difference parameter. If a single value is provided, it represents the coefficient 
and the covariate is obtained as \code{rep(c(0,1), c(ceiling(N/2), floor(N/2)))}.
}
  \item{fam}{
Which family for the response distribution? Only "gauss","pois","bin" allowed;
}
  \item{sd.eps}{
The standard deviation of the response, if \code{fam="gauss"}
}
  \item{seed}{
The \code{seed} for the random generation (used via \code{set.seed(seed)}).
}
  \item{display}{
If \code{TRUE} a (rought) plot of the simulated data (with underlying segmented profiles) is produced.
}
  \item{eta0}{
If \code{TRUE}, the true means are also returned in the dataframe.
}
} 

\details{
A list 
}
\value{
A list with components
  \item{data }{The simulated dataframe with variables, such as the response (\code{y}), the segmented covariate (\code{x}), and the `subject' id (\code{id}) if longitudinal data have been generated. The column \code{eta0} represents the true linear predictor (from which the data have been generated). The dataframe also includes the variables \code{x.diff} and \code{z.psi}, if specified. 
  }
  \item{fixef }{The fixed effects used}
  \item{psi.i }{The (potentially different) individual changepoints.}
  \item{sd}{The standard deviation of the response and of random effects, if longitudinal data have been generated}
  \item{ni }{The number of repeated measurements (per subject).}
  \item{ranef }{the random effects matrix (if longitudinal data have been generated)}
  \item{perc.psi }{the psi percentile: for each subject, the fraction of data having values smaller than the true breakpoint (if longitudinal data have been generated)}
}
\references{
Muggeo, V.M.R. (2026) Fitting segmented mixed models with random changepoints in R. \emph{J of Reproducible Statistics} \bold{1}, art 12.
10.68365/jrepstat.2026.12.
}
\author{
Vito M.R. Muggeo \email{vito.muggeo@unipa.it}
}
\note{
If \code{z.psi} and/or \code{x.diff} have been specified, no (fixed-effects) segmented relationship is drawn (even if display=TRUE) 
because it is no unique. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{segmented.lme}}
}
\examples{
\dontrun{
#simulate gaussian data: (N subjects, each with 12 measurements) 
#   following the segmented relationships  ~ a + b*x + d*(x-psi)_+
# where
# x: the segmented covariate (default is values from U(0,1))
# the intercept:        a=2 (default)
# the left slope:       b=-0.2 (default)
# the slope difference: d=0.3 (default) +.5*x.diff (x.diff built by default, see argument x.diff)
# the breakpoints is:   psi= 0.4 (default) + 0.3*Z[,1] -0.2*Z[,2]

#Each model parameter has its random effects coming from N(0, S) 

N=30
Z=cbind(1:N/N, rnorm(N))
dd=sim.seg(N,ni=12,S=diag(c(.12,.05,.3,.15)^2), 
   x.diff=.5, z.psi=list(c(.3,-.2),Z), sd.e=.003,display=TRUE)
d=dd$data

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
