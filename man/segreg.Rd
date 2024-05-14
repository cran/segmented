\name{segreg}
\alias{segreg}
\alias{stepreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting segmented regression
}
\description{
\code{segreg} (\code{stepreg}) fits linear or generalized linear segmented (stepmented) regression via a symbolic description of the linear predictor. This is an alternative, introduced since version 2.0-0 (segreg) and 2.1-0 (stepreg), to \code{segmented.lm} and \code{segmented.glm}, or (\code{stepmented.lm} or \code{stepmented.glm}).
}
\usage{
segreg(formula, data, subset, weights, na.action, family = lm, control = seg.control(), 
      transf = NULL, contrasts = NULL, model = TRUE, x = FALSE, var.psi = TRUE, ...)
      
stepreg(formula, data, subset, weights, na.action, family = lm, control = seg.control(), 
       transf = NULL, contrasts = NULL, model = TRUE, x = FALSE, var.psi = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A standard model formula also including one or more 'segmented'/'stepmented' terms via the function \code{\link{seg}}
}
  \item{data}{
The possible dataframe where the variables are stored
}
  \item{subset}{
Possible subset, as in \code{\link{lm}} or \code{\link{glm}}
}
  \item{weights}{
Possible weight vector, see \code{weights} in \code{\link{lm}} or \code{\link{glm}}
}
  \item{na.action}{
a function which indicates what happen when the data contain NA values. See \code{na.action} in \code{\link{lm}} or \code{\link{glm}}.
}
  \item{family}{
The family specification, similar to \code{family} in \code{\link{glm}}. Default to \code{'lm'} for segmented/stepmented linear models.
}
  \item{control}{
See \code{\link{seg.control}}
}
  \item{transf}{
an optional character string (with "y" as argument) meaning a function to apply to the response variable before fitting
}
  \item{contrasts}{
see \code{contrasts} in \code{\link{glm}}
}
  \item{model}{
If \code{TRUE}, the model frame is returned.
}
  \item{x}{
If \code{TRUE}, the model matrix is returned.
}
  \item{var.psi}{
logical, meaning if the standard errors for the breakpoint estimates should be returned in the object fit. If \code{FALSE}, 
the standard errors will be computed by \code{\link{vcov.segmented}} or \code{\link{summary.segmented}}. Setting \code{var.psi=FALSE} could speed up model estimation for very large datasets. Default to \code{TRUE} for \code{segreg} and \code{FALSE} for \code{stepreg}.  
}
  \item{\dots}{
Ignored
}
}
\details{
The function allows to fit segmented/stepmented (G)LM regression models using a formula interface. Results will be the same of those coming from the traditional \code{segmented.lm} and \code{segmented.glm} (or \code{stepmented.lm} or \code{stepmented.glm}), but there are some additional facilities: i) it is possible to estimate strightforwardly the segmented/stepmented relationships in each level of a categorical variable, see argument \code{by} in  \code{\link{seg}}; ii) it is possible to constrain some slopes of the segmented relationship, see argument \code{est} or \code{R} in \code{\link{seg}}.
}
\value{
An object of class "segmented" (or "stepmented") which inherits from the class "lm" or "glm" depending on \code{family} specification. See \code{\link{segmented.lm}}.
}
\references{
Muggeo, V.M.R. (2003) Estimating regression models with unknown break-points. Statistics in Medicine 22, 3055-3071.
}
\author{
Vito Muggeo
}
\note{
When the formula includes even a single segmented term with constraints (specified via the argument \code{est} in \code{seg()}), the relevant coefficients returned do not represent the slope differences as in \code{segmented.lm} or \code{segmented.glm}. The values depend on the constraints and are not usually interpretable. Use \code{\link{slope}} the recover the actual slopes of the segmented relationships. 
}

\section{Warning }{
Currently for fits returned by \code{segreg}, \code{\link{confint.segmented}} only works if \code{method="delta"}.
Constraints on the mean levels (possibly via argument 'est' of \code{seg}) are not yet allowed when calling \code{stepreg}.
} 

\seealso{
\code{\link{seg}}, \code{\link{segmented}}, \code{\link{stepmented}}
}
\examples{
###########################
#An example using segreg()
###########################

set.seed(10)
x<-1:100
z<-runif(100)
w<-runif(100,-10,-5)
y<-2+1.5*pmax(x-35,0)-1.5*pmax(x-70,0)+10*pmax(z-.5,0)+rnorm(100,0,2)

##the traditional approach
out.lm<-lm(y~x+z+w)
o<-segmented(out.lm, seg.Z=~x+z, psi=list(x=c(30,60),z=.4))

o1<-segreg(y ~ w+seg(x,npsi=2)+seg(z))
all.equal(fitted(o), fitted(o1))

#put some constraints on the slopes
o2<-segreg(y ~ w+seg(x,npsi=2, est=c(0,1,0))+seg(z))
o3<-segreg(y ~ w+seg(x,npsi=2, est=c(0,1,0))+seg(z, est=c(0,1)))

slope(o2)
slope(o3)

##see ?plant for an additional example


###########################
#An example using stepreg()
###########################

### Two stepmented covariates (with 1 and 2 breakpoints)
n=100
x<-1:n/n
z<-runif(n,2,5)
w<-rnorm(n)
mu<- 2+ 1*(x>.6)-2*(z>3)+3*(z>4)
y<- mu + rnorm(n)*.8

os <-stepreg(y~seg(x)+seg(z,2)+w) #also includes 'w' as a possible linear term
os
summary(os)
plot(os, "z", col=2:4) #plot the effect of z


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
