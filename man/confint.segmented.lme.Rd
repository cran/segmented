\name{confint.segmented.lme}
\alias{confint.segmented.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals in segmented mixed models 
}
\description{
Computes confidence intervals for all regression parameters, including the the breakpoint, in a fitted `segmented mixed' model.
}
\usage{
\method{confint}{segmented.lme}(object, parm, level = 0.95, obj.boot, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fit object returned by \code{\link{segmented.lme}}.
}
  \item{parm}{
A vector of numbers indicating which parameters should be considered. If missing all parameters.
}
  \item{level}{
The confidence level.
}
  \item{obj.boot}{
The possible list including the bootstrap distributions of the regression coefficients. Such list is returned by \code{vcov.segmented.lme(.., ret.b=TRUE)}
}
  \item{\dots}{
if \code{obj.boot} is missing and bootstrap CIs are requested, additional optional arguments, such as \code{B}, \code{seed}, and \code{it.max.b}, to be used in computations of the boot distributions. 
}
}
\details{
If \code{obj.boot} is provided or \code{...} includes the argument \code{B>0}, confidence intervals are computed by exploiting the bootstrap 
distributions.
}
\value{
A matrix (or a list of matrices if bootstrap ci are requested) including the confidence intervals for the model parameters. 
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Vito Muggeo
}
%\note{
%%  ~~further notes~~
%}

\section{Warning }{
All the functions for segmented mixed models (*.segmented.lme) are still at an experimental stage 
} 


\seealso{
 \code{\link{vcov.segmented.lme}}
}
\examples{
\dontrun{
confint(os) #asymptotic CI

confint(os, B=50) #boot CIs

#it is possible to obtain the boot distribution beforehand
ob <-vcov(os, B=50, ret.b=TRUE)
confint(os, obj.boot=ob) #boot CI

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
