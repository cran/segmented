\name{segmented}
\alias{segmented} 
\alias{segmented.default}
\alias{segmented.lm}
\alias{segmented.glm}
\alias{print.segmented}
\alias{summary.segmented}
\alias{print.summary.segmented}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Segmented relationships in regression models }
\description{
Fits regression models with segmented relationships between the response and one or more
explanatory variables. Break-point estimates are provided.
}
%- usa \cr per rompere una linea.
\usage{
segmented(obj, Z, psi, W, it.max=20, toll=0.0001, 
    visual=FALSE, last=TRUE)

%segmented.glm(obj, Z, psi, W, it.max=20, toll=0.0001, 
%    visual=FALSE, last=TRUE)
}

\arguments{
  \item{obj}{ a standard `linear' model of class \code{"lm"} or \code{"glm"}.}
  \item{Z}{ a vector or a matrix meaning the (continuous) explanatory variable(s) having segmented relationships with 
  the response of \code{obj}. }
  \item{psi}{ starting values for the break-point(s). }
  \item{W}{ an optional categorical variable specifying levels in which a segmented relationship is assumed. }
  \item{it.max}{ the maximum number of iterations allowed. }
  \item{toll}{ the tolerance value used to stop the algorithm.}
  \item{visual}{ logical value indicating whether the process fitting should be printed at each iteration. }
  \item{last}{ logical value indicating whether only the last object should be returned. }
}
\details{
Given a linear regression model (of class \code{"lm"} or \code{"glm"}), \code{segmented} 
estimates a new model having a broken-line relationship with a specified variable \code{Z}. Such a 
segmented or broken-line relationship is defined by the slope parameters and the break-points where 
the linear relation changes. Initial guesses for the break-points must be specified in the \code{psi} argument. 
The model is estimated simultaneously yielding point estimates and relevant (Wald-based) standard errors of all 
the model parameters, including the break-points.
}
\value{
The returned object depends on the \code{last} argument. If \code{last=TRUE}, the default, 
\code{segmented} returns an object of class \code{"segmented"} which inherits from the class \code{"lm"} or 
\code{"glm"} depending on the class of \code{obj}.
If \code{last=FALSE} the returned object is just a list in which the ith component is the model fitted at the ith 
approximation. Therefore the last component is the fitted model when the algorithm converges.
The function \code{summary} (i.e., \code{summary.segmented}) can be used to obtain or print a summary of the results.
  
If \code{Z} is a vector and \code{W} is missing, \code{psi} specifies the number of break-points to be estimated with 
respect to the single variable \code{Z}; if \code{Z} is a vector and 
\code{W} is provided, then a segmented relationship within each level of \code{W} is assumed: in this case 
\code{length(psi)} has to equal the number of levels of \code{W} and the ith component of \code{psi} is assumed as 
starting value of the break-point inside the ith group. 
If \code{Z} is a matrix, \code{W} has to be missing and the ith component of \code{psi} is the starting value of the 
break-point for the variable \code{Z[,i]}.

An object of class \code{"segmented"} is a list containing the components of the original object \code{obj}
with additionally the followings:
\item{psi }{Estimated break-points and relevant standard errors}
\item{it }{The number of iterations employed}
\item{epsilon }{difference in the objective function when the algorithm stops}
}
\references{ Muggeo, V.M.R. (2003) Estimating regression models with unknown break-points. 
    \emph{Statistics in Medicine} \bold{22}, 3055--3071.}
    
\author{ Vito M. R. Muggeo, \email{vito.muggeo@giustizia.it} }

\note{
\enumerate{
\item The algorithm will start if the \code{it.max} argument is greater than zero. If \code{it.max=0} 
\code{segmented} will estimate a new linear model with break-point(s) fixed at \code{psi}.
 
\item In the returned object, the name of the difference-in-slopes parameter is labelled with 
`U.name\_of\_variable'.  

\item Currently the functions can not handle missing values, thus if \code{obj} has been created with argument \code{na.action=na.omit}, \code{segmented} will not work.

\item Methods specific to the class \code{"segmented"} are 
    \itemize{
  \item \code{print.segmented}
  \item \code{summary.segmented}
  \item \code{print.summary.segmented}
            }

Others are inherited from the class \code{"lm"} or \code{"glm"} depending on the class of \code{obj}.

\item Currently multiple break-points, namely a relationships with multiple changes, are allowed only if a 
    single variable has to be analyzed, that is if \code{Z} is a vector and \code{W} is missing. 
     }
}
\section{ Warning }{It is well-known that the log-likelihood function for the break-point may be not concave, especially 
for poor clear-cut kink-relationships. In these circumstances the initial guess for the break-point, 
i.e. the \code{psi} argument, must be provided with care. For instance visual inspection of a, possibly smoothed, 
scatter-plot is usually a good way to obtain some idea on breakpoint location. 
} 

\seealso{ \code{\link{lm}}, \code{\link{glm}} }

\examples{
## A linear model with a segmented relationship in each level of 
## a categorical variable
##
## Simulate data
x<-1:100
g<-rep(0:1,c(100,100))
set.seed(15)
y1<- 2+1.5*pmax(x-60,0)+rnorm(100,0,5)
y2<- 2+0.6*pmax(x-60,0)+rnorm(100,0,5)
dati<-data.frame(xx=rep(x,2),yy=c(y1,y2),g=factor(g))
rm(x,g,y1,y2)
## Have a look at the plot
\dontrun{plot(dati$xx,dati$yy)} 
\dontrun{points(dati$xx[dati$g==0],dati$yy[dati$g==0],col=2)}
## Fit the linear model
obj<-lm(yy~0+g+xx:g,data=dati) 
## Fit segmented models
## Model I: ignore the stratification, assume an equal difference-in-slopes 
## parameter between groups
ogg0<-segmented(obj,Z=xx,psi=60)
## Model II: now stratificate by g. Here note that psi[i] refers to the 
## segmented relationship in the ith level of W
ogg1<-segmented(obj,Z=xx,W=g,psi=c(50,70),visual=TRUE)
## Results....
ogg0
summary(ogg1)
## Have a look at the fitted models
## model I
\dontrun{points(dati$xx[dati$g==0],ogg0$fitted[dati$g==0],pch=3,col=2)}
\dontrun{points(dati$xx[dati$g==1],ogg0$fitted[dati$g==1],pch=3,col=1)}
## model II
\dontrun{points(dati$xx[dati$g==0],ogg1$fitted[dati$g==0],pch=20,col=2)}
\dontrun{points(dati$xx[dati$g==1],ogg1$fitted[dati$g==1],pch=20,col=1)}
\dontrun{legend(5,60,legend=c("model I","model II"),pch=c(3,20))}
}
\keyword{ regression }
