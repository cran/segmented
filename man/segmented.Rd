\name{segmented}
\alias{segmented}
\alias{segmented.lm}
\alias{segmented.glm}
%\alias{segmented.default}
%\alias{print.segmented}
%\alias{summary.segmented}
%\alias{print.summary.segmented}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmented relationships in regression models
}
\description{
  Fits regression models with segmented relationships between the response
   and one or more explanatory variables. Break-point estimates are provided.
}
\usage{
segmented(obj, seg.Z, psi, control = seg.control(), 
    model.frame = TRUE, ...)

\method{segmented}{lm}(obj, seg.Z, psi, control = seg.control(),
    model.frame = TRUE, ...)

\method{segmented}{glm}(obj, seg.Z, psi, control = seg.control(),
    model.frame = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ standard `linear' model of class "lm" or "glm". }
  \item{seg.Z}{ a formula with no response variable, indicating the (continuous) explanatory
               variables having segmented relationships with the response
               }
  \item{psi}{ named list of vectors. The names have to match the variables of the \code{seg.Z} argument.
      Each vector includes starting values for the break-point(s) for the corresponding
      variable in \code{seg.Z} }
  \item{control}{ a list of parameters for controlling the fitting process.
      See the documentation for \code{\link{seg.control}} for details. }
  \item{model.frame}{ logical value indicating if the model.frame should be returned.}
  \item{\dots}{ optional arguments. }
}
\details{
  Given a linear regression model (of class "lm" or "glm"), segmented tries to estimate
  a new model having broken-line relationships with the variables specified in seg.Z.
  A segmented (or broken-line) relationship is defined by the slope
  parameters and the break-points where the linear relation changes. The number of breakpoints
  of each segmented relationship is fixed via the \code{psi} argument, where initial
  values for the break-points must be specified. The model
  is estimated simultaneously yielding point estimates and relevant approximate
   standard errors of all the model parameters, including the break-points.
}
\value{
  The returned object depends on the \code{last} component returned by \code{seg.control}.
  If last=TRUE, the default, segmented returns an object of class "segmented" which inherits
  from the class "lm" or "glm" depending on the class of \code{obj}. Otherwise a list is
  returned, where the last component is the fitted model at the final iteration,
  see \code{\link{seg.control}}. \cr

An object of class "segmented" is a list containing the components of the
original object \code{obj} with additionally the followings:
  \item{psi }{estimated break-points and relevant (approximate) standard errors}
  \item{it }{number of iterations employed}
  \item{epsilon }{difference in the objective function when the algorithm stops}
  \item{mframe}{ the model frame}
  \item{psi.history}{A list or a vector including the breakpoint estimates at each step}

  Other components are not of direct interest of the user.
  }

\references{ Muggeo, V.M.R. (2003) Estimating regression models with unknown 
  break-points. \emph{Statistics in Medicine} \bold{22}, 3055--3071.}
    
\author{ Vito M. R. Muggeo, \email{vmuggeo@dssm.unipa.it} }

\note{
\enumerate{
\item The algorithm will start if the \code{it.max} argument returned by \code{seg.control}
  is greater than zero. If \code{it.max=0} \code{segmented} will estimate a new linear model with
 break-point(s) fixed at the values reported in \code{psi}.
 
\item In the returned object, the name of the difference-in-slopes parameter is 
labelled with `U.name\_of\_variable'.

\item Methods specific to the class \code{"segmented"} are
    \itemize{
  \item \code{print.segmented}
  \item \code{summary.segmented}
  \item \code{print.summary.segmented}
            }

Others are inherited from the class \code{"lm"} or \code{"glm"} depending on the
 class of \code{obj}.

     }
}
\section{ Warning }{It is well-known that the log-likelihood function for the 
break-point may be not concave, especially 
for poor clear-cut kink-relationships. In these circumstances the initial guess
 for the break-point, 
i.e. the \code{psi} argument, must be provided with care. For instance visual 
inspection of a, possibly smoothed, 
scatter-plot is usually a good way to obtain some idea on breakpoint location. 
Moreover it is also advisable to look at the coefficients of the `gap' variables. At the
convergence they should be small or at least, non-significantly different from zero.
\code{\link{summary.segmented}} and \code{\link{print.summary.segmented}} return
information of the `gap' coefficients.
}

\seealso{ \code{\link{lm}}, \code{\link{glm}} }

\examples{

set.seed(12)
xx<-1:100
zz<-runif(100)
yy<-2+1.5*pmax(xx-35,0)-1.5*pmax(xx-70,0)+15*pmax(zz-.5,0)+rnorm(100,0,2)
dati<-data.frame(x=xx,y=yy,z=zz)
out.lm<-lm(y~x,data=dati)
o<-\method{segmented}{lm}(out.lm,seg.Z=~x,psi=list(x=c(30,60)),
    control=seg.control(display=FALSE))
slope(o)

out.lm1<-lm(y~z,data=dati)
o1<-update(o,seg.Z=~x+z,psi=list(x=c(30,60),z=.3))

            }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{nonlinear }
