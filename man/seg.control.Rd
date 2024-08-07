\name{seg.control}
\alias{seg.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary for controlling segmented/stepmented model fitting }
\description{
  Auxiliary function as user interface for 'segmented' and 'stepmented' fitting. Typically
     only used when calling any 'segmented' or 'stepmented' method.
}
\usage{
seg.control(n.boot=10, display = FALSE, tol = 1e-05, it.max = 30, fix.npsi=TRUE,
    K = 10, quant = FALSE, maxit.glm = NULL, h = 1.25, break.boot=5, size.boot=NULL, 
    jt=FALSE, nonParam=TRUE, random=TRUE, seed=NULL, fn.obj=NULL, digits=NULL, 
    alpha = NULL, fc=.95, check.next=TRUE, tol.opt=NULL, fit.psi0=NULL, eta=NULL, 
    min.nj=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.boot}{ number of bootstrap samples used in the bootstrap restarting algorithm. If 0 the standard algorithm,
    i.e. without bootstrap restart, is used. Default to 10 that appears to be sufficient in most of problems. However
    when multiple breakpoints have to be estimated it is suggested to increase \code{n.boot}, e.g. \code{n.boot=50}, and even \code{break.boot}.}
  \item{display}{ logical indicating if the value of the objective function should be printed along with current breakpoint estimates at each iteration or at each bootstrap resample (but no more than 5 breakpoints are printed). If bootstrap restarting is employed, the values of objective and breakpoint estimates should not change at the last runs.}
  \item{tol}{ positive convergence tolerance. }
  \item{it.max}{ integer giving the maximal number of iterations. }
  \item{fix.npsi}{logical (it replaces previous argument \code{stop.if.error}) If \code{TRUE} (default) the \emph{number} (and not location) of breakpoints is held fixed           throughout iterations. Otherwise a sort of `automatic' breakpoint selection is carried out, provided that several starting values are supplied for the breakpoints, 
    see argument \code{psi} in \code{\link{segmented.lm}} or \code{\link{segmented.glm}}. The idea, relying on removing the `non-admissible' breakpoint estimates at each     iteration, is discussed in Muggeo and Adelfio (2011) and it is not compatible with the bootstrap restart algorithm. \code{fix.npsi=FALSE}, indeed, should be considered as a preliminary and tentative approach to deal with an unknown number of breakpoints. Alternatively, see \code{\link{selgmented}}.} 
  \item{K}{ the number of quantiles (or equally-spaced values) to supply as starting values for the breakpoints 
    when the \code{psi} argument of \code{segmented} is set to \code{NA}. \code{K} is ignored when \code{psi} is different from \code{NA}.
    }
  \item{quant}{logical, indicating how the starting values should be selected. If \code{FALSE} equally-spaced 
  values are used, otherwise the quantiles. Ignored when \code{psi} is different from \code{NA}.} 
  \item{maxit.glm}{ integer giving the maximum number of inner IWLS iterations (see details). If \code{NULL}, the number is low in the first iterations and then increases as th eprocess goes on. Ignored for segmented lm fits }
  \item{h}{ positive factor modifying the increments in breakpoint updates during
    the estimation process (see details). }
  \item{break.boot}{
  Integer, less than \code{n.boot}. If \code{break.boot} consecutive bootstrap samples lead to the same objective function during the estimation process, the algorithm stops without performing all \code{n.boot} 'trials'.
  This can save computational time considerably. Default is \code{5} for the \code{segmented} and \code{5+2} for the \code{stepmented} functions. However if the number of changepoints is large, \code{break.boot} should be increased, even 10 or 15.
    }
  \item{size.boot}{the size of the bootstrap samples. If \code{NULL}, it is taken equal to the actual sample size. If the sample is very large, the idea is to run bootstrap restarting using smaller bootstrap samples.}
  \item{jt}{logical. If \code{TRUE} the values of the segmented variable(s) are jittered before fitting the model to the
  bootstrap resamples.}
  \item{nonParam}{
  if \code{TRUE} nonparametric bootstrap (i.e. case-resampling) is used, otherwise residual-based.
  Currently working only for LM fits. It is not clear what residuals should be used for GLMs.}
  \item{random}{
  if \code{TRUE}, when the algorithm fails to obtain a solution, random values are employed to obtain candidate values.
  }
  \item{seed}{
  The seed to be passed on to \code{set.seed()} when \code{n.boot>0}. If \code{NULL}, a seed depending on the response values is generated and used. Otherwise it can be a numerical value or, if \code{NA}, a random value is generated. 
  Fixing the seed can be useful to replicate \emph{exactly} the results when the bootstrap restart algorithm is employed. Whichever choice, the segmented fit includes the component \code{seed} representing the value saved just before the bootstrap resampling. Re-use it if you want to replicate the bootstrap restarting algorithm with the \emph{same} re-samples.
  }
  \item{fn.obj}{
  A \emph{character string} to be used (optionally) only when \code{segmented.default} is used. It represents the function 
  (with argument \code{'x'}) to be applied to the fit object to extract the objective function to be \emph{minimized}. 
  Thus for \code{"lm"} fits (although unnecessary) it should be \code{fn.obj="sum(x$residuals^2)"}, for
  \code{"coxph"} fits it should be \code{fn.obj="-x$loglik[2]"}. If \code{NULL} the `minus log likelihood' extracted from 
  the object, namely \code{"-logLik(x)"}, is used. See \code{\link{segmented.default}}.
    }
  \item{digits}{optional. If specified it means the desidered number of decimal points of the breakpoint to be used during the iterative algorithm.
  }
%  \item{conv.psi}{optional. Should convergence of iterative procedure to be assessed on changes of breakpoint estimates or changes in the %objective? Default to FALSE.
%  }
  \item{alpha}{optional numerical values. The breakpoints are estimated within the quantiles \code{alpha[1]} and \code{alpha[2]} of the relevant covariate. If a single value is provided, it is assumed \code{alpha} and \code{1-alpha}. Defaults to \code{NULL} which means \code{alpha=max(.05, 1/n)}. Note: Providing \code{alpha=c(mean(x<=a),mean(x<=b))} means to constrain the breakpoint estimates within \eqn{[a,b]}{[a,b]}.
  }
%  \item{min.step}{optional. The minimum step size to break the iterative algorithm. Default to 0.0001.}
%  \item{stop.if.error}{ same than \code{fix.npsi}. \emph{This argument will be removed in next releases}, %and replaced by 
%  \code{fix.npsi}.
% If provided, and different from \code{NULL}, it overwrites \code{fix.npsi}
%  }
  \item{fc}{A proportionality factor (\eqn{\le 1}{<= 1}) to adjust the breakpoint estimates \emph{if} these come close to the boundary or too close each other. For instance, if \code{psi} turns up close to the maximum, it will be changed to \code{psi*fc} or to \code{psi/fc} if close to the minimum. This is useful to get finite point estimate and standard errors for each slope paramete.
  }
  \item{check.next}{logical, effective only for stepmented fit. If \code{TRUE} the solutions next to the current one are also investigated.
  }
  \item{tol.opt}{Numerical value to be passed to \code{tol} in \code{\link{optimize}}.
  }
  \item{fit.psi0}{Possible list including preliminary values.
  }
  \item{eta}{Only for segmented/stepmented fits: starting values to be passed to \code{etastart} in \code{\link{glm.fit}}.
  }
  \item{min.nj}{How many observations (at least) should be in the covariate intervals induced by the breakpoints?}
  }
\details{
  Fitting a `segmented' GLM model is attained via fitting iteratively standard GLMs. The number of (outer) 
  iterations is governed by \code{it.max}, while the (maximum) number of (inner) iterations to fit the GLM at 
  each fixed value of psi is fixed via \code{maxit.glm}. Usually three-four inner iterations may be sufficient.

  When the starting value for the breakpoints is set to \code{NA} for any segmented variable specified
  in \code{seg.Z}, \code{K} values (quantiles or equally-spaced) are selected as starting values for the breakpoints. 
  %In this case, it may be useful to set also \code{fix.npsi=FALSE} to automate the procedure, see Muggeo and Adelfio (2011).
  %The maximum number of iterations (\code{it.max}) should be also increased when the `automatic' procedure is used.

%  If \code{last=TRUE}, the object resulting from \code{segmented.lm} (or \code{segmented.glm}) is a
%  list of fitted GLM; the i-th model is the segmented model with the values of the breakpoints at the i-th %iteration. 

  Since version 0.2-9.0 \code{segmented} implements the bootstrap restarting algorithm described in Wood (2001).
  The bootstrap restarting is expected to escape the local optima of the objective function when the
  segmented relationship is noisy and the loglikelihood can be flat. Notice bootstrap restart runs \code{n.boot} iterations regardless of \code{tol} that only affects convergence within the inner loop.
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
     #seg.control(display = TRUE, maxit.glm=4)
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
