\name{summary.segmented.lme}
\alias{summary.segmented.lme}
%\alias{print.summary.segmented}
\title{ Summarizing  model fits for segmented regression }
\description{
  summary method for class \code{segmented.lme}.
}
\usage{
\method{summary}{segmented.lme}(object, .vcov=NULL, digits = max(3, getOption("digits") - 3), ...)

%\method{print}{summary.segmented}(x, short=x$short, var.diff=x$var.diff, 
%    digits = max(3, getOption("digits") - 3),
%    signif.stars = getOption("show.signif.stars"),...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class "segmented.lme". }
  
  \item{.vcov}{ Optional. The full covariance matrix for the parameter estimates. If provided, standard errors are computed (and displayed) according to this matrix.}
%  \item{x}{a \code{summary.segmented} object produced by \code{summary.segmented()}.}
  \item{digits}{controls number of digits printed in output.}
%  \item{signif.stars}{logical, should stars be printed on summary tables of coefficients?}
  \item{\dots}{ further arguments. }
}
\details{
  %If \code{short=TRUE} only coefficients of the segmented relationships are printed.
  Da finire
  
}
\value{
A list (similar to one returned by \code{segmented.lm}) with additional components:
  \item{psi}{estimated break-points and relevant (approximate) standard errors}
%  \item{Ttable }{estimates and standard errors of the model parameters. This is similar
%  to the matrix \code{coefficients} returned by \code{summary.lm} or \code{summary.glm},
%  but without the rows corresponding to the breakpoints. Even the p-values relevant to the
%  difference-in-slope parameters have been replaced by NA, since they are meaningless in
%  this case, see \code{\link{davies.test}}.}
%  \item{cov.var.diff}{if \code{var.diff=TRUE}, the covaraince matrix accounting for heteroscedastic errors.}
%  \item{sigma.new}{if \code{var.diff=TRUE}, the square root of the estimated error variances in each interval.}
%  \item{df.new}{if \code{var.diff=TRUE}, the residual degrees of freedom in each interval.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Vito M.R. Muggeo }
\seealso{ \code{\link{print.segmented}} }
\examples{
##continues example from segmented.lme()
# summary(os)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
