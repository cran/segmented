\name{segmented.lme}
\alias{segmented.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmented relationships in linear mixed models
}
\description{
Fits linear mixed models with a segmented relationship between the response and 
a numeric covariate. Random effects are allowed in each model parameter, including the breakpoint.
}
\usage{
\method{segmented}{lme}(obj, seg.Z, psi, npsi = 1, fixed.psi = NULL, 
    control = seg.control(), model = TRUE, 
    z.psi = ~1, x.diff = ~1, random = NULL, 
    random.noG = NULL, start.pd = NULL, psi.link = c("identity", "logit"), 
    nq = 0, adjust = 0, start = NULL, data, fixed.parms = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A 'lme' fit returned by \code{lme} or simply its call. See example below. This represents the linear mixed model
where the segmented relationship is added.
}
  \item{seg.Z}{
A one-sided formula indicating the segmented variable, i.e. the quantitative variable having a segmented relationship with
  the response. In longitudinal studies typically it is the time.}
  \item{psi}{
An optional starting value for the breakpoint. If missing a starting value is obtained via the nadir estimate of a quadratic fit.
When provided it may be a single numeric value or a vector of length equal to the number of clusters (i.e. subjects).
}
  \item{z.psi}{
Optional. A one-sided formula meaning the sub-regression model for the changepoint parameter. Default to \code{~1}. 
}
  \item{x.diff}{
Optional. A one-sided formula meaning the sub-regression model for the difference-in-slopes parameter. 
Default to \code{~1} for no covariate for the difference-in-slope parameter.
}
  \item{npsi}{
Ignored. Currently only \code{npsi=1} is allowed.
}
  \item{fixed.psi}{
Ignored.
}
  \item{control}{
A list returned by \code{\link{seg.control}}, in particular \code{display}, \code{n.boot} for the bootstrap restarting.
}
  \item{model}{
Ignored.
}
  \item{random}{
A list, as the one supplied in \code{random} of \code{lme()} including the random effects. Default to \code{NULL},
meaning that the same random effect structure of the initial lme fit supplied in \code{obj} should be used. 
When specified, this list could include the variables '\code{G0}' and '\code{U}'. 
\code{G0} means random effects in the breakpoints and \code{U} means random effects in the
slope-difference parameter. Assuming \code{id} is the the cluster variable and \code{x} the segmented variable, some examples are


\code{random = list(id = pdDiag(~1 + x + U))} #ind. random eff. (changepoint fixed)

\code{random = list(id = pdDiag(~1 + x + U + G0))} #ind. random eff. (in the changepoint too) 

\code{random = list(id=pdBlocked(list(pdSymm(~1+x), pdSymm(~U+G0-1))))} #block diagonal 
}
  \item{random.noG}{
Ignored.
}
  \item{start.pd}{
An optional starting value for the variances of the random effects. It should be coherent with the
  specification in \code{random}.
}
  \item{psi.link}{
The link function used to specify the sub-regression model for the breakpoint \eqn{psi}.
The identity (default) assumes
\deqn{\psi_i=\eta_i}
while the logit link is
\deqn{\psi_i=(m+M*exp(\eta_i))/(1+exp(\eta_i))}
where m and M are the observed minimum and maximum of the segmented variable \code{Z}.
In each case the `linear predictor' is \eqn{\eta_i=\kappa_0+\kappa_1*z_i+k_i}, where z is 
the covariate specified in \code{z.psi} and the \eqn{k_i}s are the changepoint random effects 
included by means of \code{G0} in the \code{random} argument.
}
  \item{nq}{
Integer value to assess if the estimated breakpoint for the subject i can be considered 'reliable' or not. For each subject, \eqn{\hat\psi} is not considered reliable if 
\deqn{\hat{\psi}\le z_{(1+nq)} \; \mathrm{ and } \; \hat{\psi}\ge z_{(n_i-nq)} }
where \eqn{z_{(\cdot)}} means the sorted values of the segmented variable. Note \code{nq} just affects the logical vector \code{attr(psi.i,"is.break")} where \code{psi.i} is the component of the model fit including the subject-specific breakpoints. Default is \code{nq=0}.
}
  \item{adjust}{
A numerical (0 or 1) value. If \code{adjust=1} at convergence and for each subject \code{i} the estimated breakpoint assessed as unreliable (see argument \code{nq}) is moved to the maximum of the covariate range (for the subject i). Thus for the subjects with 'unreliable' breakpoint, the corresponding fitted  profile will be truly linear. 
}
  \item{start}{
An optional list including the \emph{starting values} for the difference-in-slopes parameter, delta0 and delta,
and the changepoint parameter, kappa and kappa0. When provided, 'kappa0' overwrites 'psi'.

If provided, the components 'delta' and 'kappa' should be \emph{named} vectors with length and names
matching length and names in \code{x.diff} and \code{z.psi} respectively. The component \code{delta0} can be
a scalar or a vector with length equal to the number of clusters (subjects).
}
  \item{data}{
the dataframe where the variables are stored. If missing, the dataframe of the \code{"lme"} fit \code{obj} is assumed.
}
  \item{fixed.parms}{
An optional \emph{named} vector representing the coefficients \emph{of the changepoint} to be maintained \emph{fixed}
during the estimation process. Allowed names are "G0" or any variable (in the dataframe) supposed to
affect the location of breakpoints.
For instance \code{fixed.parms=c(G0=.3)} implies a fixed value for the changepoint.
Notice if you use the same variable in \code{fixed.parms} and in \code{z.psi}, for instance
\code{fixed.parms=c(x2=.3)} and \code{z.psi=~x2}, a warning is printed and the
coefficient "G.x2" is estimated to maximize the log likelihood \emph{given} that fixed value. As an example, suppose the 
unconstrained estimated coefficient for x2, say, in \code{z.psi} is 0.5; if in a new call both \code{fixed.parms=c(x2=.4)} and 
\code{z.psi=~x2} are included, the estimate of "G.x2" will be (approximately) 0.1.
Essentially, if you really want to fix the parameters in \code{fixed.parms}, then 
do not include the same covariates in \code{z.psi}.
}
  \item{...}{
Ignored
}

%  \item{tol}{
%the tolerance value to declare convergence
%}
%  \item{it.max}{
%the maximum number of iterations allowed
%}
%  \item{display}{
%logical. If \code{TRUE} the objective function value is printed at each iteration
%}
}
\details{
The function fits segmented mixed regression models, i.e. segmented models with random effects also in the slope-difference and change-point parameters. 
}
\value{
A list of class \code{segmented.lme} with several components. The most relevant are 
  \item{lme.fit }{The fitted lme object at convergence}
  \item{lme.fit.noG }{The fitted lme object at convergence assuming known the breakpoints}
  \item{psi.i}{The subject/cluster-specific change points (fixed + random). It includes 2 attributes: \code{attr(,"ni")} for the number of measurements in each 'cluster', and \code{attr(,"is.break")} a vector of logicals indicating if the breakpoint for each subject i can be reliable (\code{TRUE}) or not (\code{FALSE}). Here 'reliable' simply means within the covariate range (for subject i). See also argument \code{nq}.}
  \item{fixed.eta.psi}{The fixed-effect linear predictor for the change points regression equation. These values will different among 'clusters' only if at least one covariate has been specified in \code{z.psi}.}
  \item{fixed.eta.delta}{The fixed-effect linear predictor of the slope difference regression equation. These values will different among 'clusters' only if at least one covariate has been specified in \code{x.diff}.}
}
\references{
Muggeo V., Atkins D.C., Gallop R.J., Dimidjian S. (2014) 
Segmented mixed models with random changepoints: a maximum likelihood approach 
with application to treatment for depression study. Statistical Modelling, 14, 293-313.

Muggeo V. (2016) 
Segmented mixed models with random changepoints in R.  
Working paper available on RG. doi: 10.13140/RG.2.1.4180.8402
}
\author{
Vito M.R. Muggeo \email{vito.muggeo@unipa.it}
}
\note{
Currently only one breakpoint (with or without random effects) can be estimated.
}

\section{Warning }{
All the functions for segmented mixed models (*.segmented.lme) are still at an experimental stage 
} 


\seealso{
\code{\link{plot.segmented.lme}} for the plotting method and \code{\link{segmented.default}} (example 2) for segmented models with no random effects in breakpoints or slope difference.
}
\examples{

\dontrun{
library(nlme)
data(Cefamandole)
Cefamandole$lTime <-log(Cefamandole$Time)
Cefamandole$lconc <-log(Cefamandole$conc)

o<-lme(lconc ~ lTime, random=~1|Subject, data=Cefamandole)

os<-segmented.lme(o, ~lTime, random=list(Subject=pdDiag(~1+lTime+U+G0)), 
  control=seg.control(n.boot=0, display=TRUE))
slope(os)


####################################################
# covariate effect on the changepoint and slope diff


#let's assume a new subject-specific covariates..
set.seed(69)
Cefamandole$z <- rep(runif(6), rep(14,6))
Cefamandole$group <- gl(2,42,labels=c('a','b'))

#Here 'group' affects the slopes and 'z' affects the changepoint 

o1 <-lme(lconc ~ lTime*group, random=~1|Subject, data=Cefamandole)
os1 <- segmented(o1, ~lTime, x.diff=~group, z.psi=~z, 
  random=list(Subject=pdDiag(~1+lTime+U+G0)))

slope(os1, by=list(z=.1)) #the slope estimates when z=.1


###################################################
# A somewhat "complicated" example:
#     i)  strong heterogeneity in the changepoints
#     ii) No changepoint for the Subject #7 (added) 

d<-Cefamandole
d$x<- d$lTime
d$x[d$Subject==1]<- d$lTime[d$Subject==1]+3
d$x[d$Subject==5]<- d$lTime[d$Subject==5]+5
d$x[d$Subject==3]<- d$lTime[d$Subject==3]-5
d<-rbind(d, d[71:76,])
d$Subject <- factor(d$Subject, levels=c(levels(d$Subject),"7")) 
d$Subject[85:90] <- rep("7",6)

o<-lme(lconc ~ x, random=~1|Subject, data=d)
os2<-segmented.lme(o, ~x, random=list(Subject=pdDiag(~1+x+U+G0)), 
  control=seg.control(n.boot=5, display=TRUE))

plot(os2, n.plot = c(3,3))
os2$psi.i
attr(os2$psi.i, "is.break") #it is FALSE for Subject #7

plot(os2, n.plot = c(3,3), xscale=-1, yscale = -1) #plots with subject-specific scales
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
