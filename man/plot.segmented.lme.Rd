\name{plot.segmented.lme}
\alias{plot.segmented.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for segmented mixed objects
}
\description{
Takes a fitted \code{segmented.lme} object returned by \code{segmented()} and plots (or adds)
  the fitted broken-line relationship for the selected segmented term.
}
\usage{
\method{plot}{segmented.lme}(x, level=1, id = NULL, res = TRUE, pop = FALSE, yscale = 1, xscale = 1, 
    n.plot, pos.leg = "topright", vline = FALSE, lines = TRUE, 
    by=NULL, add=FALSE, conf.level=0, withI=TRUE, vcov.=NULL, shade=FALSE, 
    drop.var=NULL, text.leg=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{"segmented.lme"}
}
  \item{level}{
An integer giving the level of grouping to be used when computing the segmented relationship(s). \code{level=0} means depending on fixed effects estimates only (such estimates are also said, to some extend, 'population' or 'marginal' estimates), otherewise the segmented lines will also depend on the random effects predictions. 
}

  \item{id}{
A scalar or vector meaning which subjects profiles have to be drawn. If \code{NULL} (default) all profiles are drawn. Ignored if \code{level=0}.
}
  \item{res}{
If \code{TRUE} the data points are also drawn. Ignored if \code{level=0}.
}
  \item{pop}{
if \code{TRUE}, the fitted segmented relationships based on fixed-effects only is also portrayed. Ignored if \code{level=0}.
}
  \item{yscale}{
If \code{>= 0}, the same and common y-scale is used for all 'subjects' (panels); otherwise the y-scale will depend on the actual (observed and fitted) values for each 'subject'.
}
  \item{xscale}{
If \code{>= 0}, the same and common x-scale is used for all 'subjects' (panels); otherwise the x-scale will depend on the actual observed values of the segmented covariate for each 'subject'.
}
  \item{n.plot}{
a vector to be passed to \code{par(mfrow=c(..))} for plotting multiple panels (should be coherent with \code{length(id)}). If missing, it is computed automatically depending on \code{length(id)}. 
}
  \item{pos.leg}{
a character ('topright', 'topleft',...) meaning the location of the legend. Set \code{NULL} for no legend.  
}
  \item{vline}{
logical, if \code{TRUE} a vertical dashed segment is added to emphasize the breakpoint location.
}
  \item{lines}{
logical, if \code{FALSE} points, rather than lines, are used to portray the segmented relationships.
}
  \item{by}{
A named list indicating covariate names and corresponding values affecting the fitted segmented relationship. For instance: 
\code{by=list(sex="male",z=.2)}, provided that the variables \code{sex} and \code{z} affect the segmented relationship. Effective only if \code{level=0}.
}
  \item{add}{
If \code{TRUE} the (fixed-effect) fitted segmented relationship is added to the current plot. Effective only if \code{level=0}. 
}
  \item{conf.level}{
The confidence level for pointwise confidence intervals. Effective only if \code{level=0}.
}
  \item{withI}{
If \code{TRUE}, the level 0 segmented relationship is computed with the model intercept. Effective only if \code{level=0}.
}
  \item{vcov.}{
The fixed effects covariance matrix. If \code{NULL}, it is computed by \code{vcov.segmented.lme()}. Effective only if \code{level=0}.
}
  \item{shade}{
If \code{TRUE}  (and \code{conf.level>0}) the area within the pointiwise CIs is shaded. Effective only if \code{level=0}.
}
  \item{drop.var}{
Possible coefficient names to be removed before computing the segmented relationship (E.g. the group-specific intercept.).
}
  \item{text.leg}{
If specified (and \code{pos.leg} has been also specified), it is the legend text. Effective only if \code{level=0}.
}
  \item{\dots}{
additional arguments, such as 'ylab', 'xlab', 'col.l' (for the fitted line) and 'col' (for the data points, provided that \code{res=TRUE} has been set). If \code{level=0} and \code{conf.level>0}, \code{lty} and \code{lwd} can be vectors.
}
}
\details{
The function plots the 'subject'-specific segmented profiles for the 'subjects' specificed in \code{id} or, if \code{level=0}, the fitted segmented relationship based on fixed effects only.
}
\value{
A single or multiple (depending on \code{level} and \code{id}) plot showing the fitted segmented profiles.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
\note{
  If \code{by} is specified (and \code{level=0} is set), a legend is also added in the plot reporting covariate(s) name and value affecting the segmented relationship.  Set \code{pos.leg=TRUE} to have no legend. On the other hand, use \code{text.leg} to add legend reporting the covariate baseline values.
}

\section{Warning }{
All the functions for segmented mixed models (*.segmented.lme) are still at an experimental stage 
} 


\seealso{
\code{\link{segmented.lme}}
}
\examples{
\dontrun{
#continues example from segmented.lme

plot(os, yscale=-1) #different y-scales

}
}
\keyword{ regression }
\keyword{ nonlinear }
\concept{changepoint}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
