\name{seg.control}
\alias{seg.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary for controlling segmented model fitting }
\description{
  Auxiliary function as user interface for 'segmented' fitting. Typically
     only used when calling any 'segmented' method (\code{segmented.lm} or \code{segmented.glm}).
}
\usage{
seg.control(toll = 1e-04, it.max = 10, display = FALSE,
    stop.if.error = TRUE, K = 10, quant = FALSE, last = TRUE, maxit.glm = 25, h = 1, 
    n.boot=20, size.boot=NULL, gap=FALSE, jt=FALSE, nonParam=TRUE,
    random=TRUE, powers=c(1,1), seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{toll}{ positive convergence tolerance. }
  \item{it.max}{ integer giving the maximal number of iterations. }
  \item{display}{ logical indicating if the value of the \emph{working} objective function should be printed
    at each iteration. The \emph{working} objective function is the objective function of the working model
    including the gap coefficients.}
  \item{stop.if.error}{ logical indicating if non-admissible break-points should be removed 
    during the estimating algorithm. Set it to \code{FALSE} if you want to perform a sort of
    `automatic' breakpoint selection, provided that several starting values are provided for the breakpoints.
    See argument \code{psi} in \code{\link{segmented.lm}} or \code{\link{segmented.glm}}. 
    The idea of removing `non-admissible' break-points during the iterative process is discussed in 
    Muggeo and Adelfio (2011) and it is not compatible with the bootstrap restart algorithm.
    This approach, indeed, should be considered as a preliminary and tentative approach to deal with an unknown number of 
    breakpoints.
    }
  \item{K}{ the number of quantiles (or equally-spaced values) to supply as starting values for the breakpoints 
    when the \code{psi} argument of \code{segmented} is set to \code{NA}. 
    \code{K} is ignored when \code{psi} is different from \code{NA}.
    }
  \item{quant}{logical, indicating how the starting values should be selected. If \code{FALSE} equally-spaced 
  values are used, otherwise the quantiles. Ignored when \code{psi} is different from \code{NA}.} 
  \item{last}{ logical indicating if output should include only the last fitted model.}
  \item{maxit.glm}{ integer giving the maximum number of inner IWLS iterations (see details). }
  \item{h}{ positive factor (from zero to one) modifying the increments in breakpoint updates during
    the estimation process (see details). }
  \item{n.boot}{ number of bootstrap samples used in the bootstrap restarting algorithm. If 0 the standard algorithm,
    i.e. without bootstrap restart, is used. Default to 10 that appears to be sufficient in most of problems. However
    when multiple breakpoints have to be estimated it is suggested to increase \code{n.boot}, e.g. \code{n.boot=50}.}
  \item{size.boot}{the size of the bootstrap samples. If \code{NULL}, it is taken equal to the actual sample size.}
  \item{gap}{logical, if \code{FALSE} the gap coefficients are \emph{always} constrained to zero at the convergence.}
  \item{jt}{logical. If \code{TRUE} the values of the segmented variable(s) are jittered before fitting the model to the
  bootstrap resamples.}
  \item{nonParam}{
  if \code{TRUE} nonparametric bootstrap (i.e. case-resampling) is used, otherwise residual-based.
  Currently working only for LM fits. It is not clear what residuals should be used for GLMs.}
  \item{random}{
  if \code{TRUE}, when the algorithm fails to obtain a solution, random values are employed to obtain candidate values.
  }
  \item{powers}{
  The powers of the pseudo covariates employed by the algorithm. These are possibly altered during the iterative 
  process to stabilize the estimation procedure. Usually of no interest for the user.
  }
  \item{seed}{
  The seed to be passed on \code{set.seed()} when \code{n.boot>0}. Setting the seed can be useful to replicate
  the results when the bootstrap restart algorithm is employed. In fact \code{employed.Random.seed} includes 
  the integer vector saved just before the bootstrap resampling. Use it if you want to replicate the bootstrap
  restarting algorithm with the \emph{same} samples.
  }
  }

\details{
  Fitting a `segmented' GLM model is attained via fitting iteratively standard GLMs. The number of (outer) 
  iterations is governed by \code{it.max}, while the (maximum) number of (inner) iterations to fit the GLM at 
  each fixed value of psi is fixed via \code{maxit.glm}. Usually three-four inner iterations may be sufficient.

  When the starting value for the breakpoints is set to \code{NA} for any segmented variable specified
  in \code{seg.Z}, \code{K} values (quantiles or equally-spaced) are selected as starting values for the breakpoints. 
  In this case, it may be useful to set also \code{stop.if.error=FALSE} to automate the procedure, see Muggeo and Adelfio (2011).
  The maximum number of iterations (\code{it.max}) should be also increased when the `automatic' procedure is used.

  If \code{last=TRUE}, the object resulting from \code{segmented.lm} (or \code{segmented.glm}) is a
  list of fitted GLM; the i-th model is the segmented model with the values of the breakpoints at the i-th iteration. 

  Sometimes to stabilize the procedure, it can be useful to set \code{h<1} to reduce the increments in the breakpoint
  updates. At each
  iteration the updated estimate is usually given by \code{psi.new=psi.old+increm}. By setting \code{h<1}
  (actually \code{min(abs(h),1)} is considered) causes the following updates of the breakpoint estimate:
  \code{psi.new=psi.old+h*increm}.

  Since version 0.2-9.0 \code{segmented} implements the bootstrap restarting algorithm described in Wood (2011).
  The bootstrap restarting is expected to escape the local optima of the objective function when the
  segmented relationship is flat. Notice bootstrap restart runs \code{n.boot} iterations regardless of 
  \code{toll} that only affects convergence within the inner loop.
}
\value{
  A list with the arguments as components.
}
\references{ Muggeo, V.M.R., Adelfio, G. (2011) Efficient change point detection in genomic sequences of 
    continuous measurements. \emph{Bioinformatics} \bold{27}, 161--166. 

Wood, S. N. (2001) Minimizing model fitting objectives that contain spurious local minima
    by bootstrap restarting. \emph{Biometrics} \bold{57}, 240--244. 
    }

\author{ Vito Muggeo }
%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

     #decrease the maximum number inner iterations and display the 
     #evolution of the (outer) iterations
     seg.control(display = TRUE, maxit.glm=4)
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
