\name{summary.segmented}
\alias{summary.segmented}
\alias{print.summary.segmented}
\title{ Summarizing  model fits for segmented regression }
\description{
  summary method for class \code{segmented}.
}
\usage{
\method{summary}{segmented}(object, short = FALSE, var.diff = FALSE, p.df="p", .vcov=NULL, ...)

\method{print}{summary.segmented}(x, short=x$short, var.diff=x$var.diff, 
    digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"),...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class "segmented". }
  \item{short}{ logical indicating if the `short' summary should be printed. }
  \item{var.diff}{ logical indicating if different error variances should be computed 
    in each interval of the segmented variable, see Details. If \code{.vcov} is provided, \code{var.diff} is set to \code{FALSE}. }
  \item{p.df}{ A character as a function of \code{'p'} (number of parameters) and \code{'K'} (number of groups or segments) affecting computations of the group-specific 
  variance (and the standard errors) if \code{var.diff=TRUE}, see Details.}
  \item{.vcov}{ Optional. The full covariance matrix for the parameter estimates. If provided, standard errors are computed (and displayed) according to this matrix.}
  \item{x}{a \code{summary.segmented} object produced by \code{summary.segmented()}.}
  \item{digits}{controls number of digits printed in output.}
  \item{signif.stars}{logical, should stars be printed on summary tables of coefficients?}
  \item{\dots}{ further arguments. }
}
\details{
  If \code{short=TRUE} only coefficients of the segmented relationships are printed.
  If \code{var.diff=TRUE} and there is only one segmented variable, different error variances are 
  computed in the intervals defined by the estimated breakpoints of the segmented variable. 
  For the jth interval with \eqn{n_j}{nj} observations, the error variance is estimated via \eqn{RSS_j/(n_j-p)}{RSSj/(nj-p)}, 
  where \eqn{RSS_j} is the residual sum of squares in interval j, and \eqn{p}{p} is the number of model parameters. This number to be subtracted from \eqn{n_j}{nj} can be changed via argument \code{p.df}. For instance \code{p.df="0"} uses \eqn{RSS_j/(n_j)}{RSSj/(nj)}, and \code{p.df="p/K"} leads to \eqn{RSS_j/(n_j-p/K)}{RSSj/(nj-p/K)}, where \eqn{K}{K} is the number of groups (segments), and \eqn{p/K}{p/K} can be interpreted as the average number of model parameter in that group. 
  
  Note \code{var.diff=TRUE} only affects the estimates covariance matrix. It does \emph{not} affect the parameter estimates, neither the log likelihood and relevant measures, such as AIC  or BIC. In other words, \code{var.diff=TRUE} just provides 'alternative' standard errors, probably appropriate when the error variances are different before/after the estimated breakpoints. Also \eqn{p-values}{p-values} are computed using the t-distribution with 'naive' degrees of freedom (as reported in \code{object$df.residual}).
  
  If \code{var.diff=TRUE} the variance-covariance matrix of the estimates is computed via the 
  sandwich formula, 
  \deqn{(X^TX)^{-1}X^TVX(X^TX)^{-1}}{(X'X)^{-1}X'VX(X'X)^{-1}} 
  where V is the diagonal matrix including the different group-specific error variance estimates. Standard errors are the square root of the main diagonal of this matrix.
}
\value{
A list (similar to one returned by \code{segmented.lm} or \code{segmented.glm}) with additional components:
  \item{psi }{estimated break-points and relevant (approximate) standard errors}
  \item{Ttable }{estimates and standard errors of the model parameters. This is similar
  to the matrix \code{coefficients} returned by \code{summary.lm} or \code{summary.glm},
  but without the rows corresponding to the breakpoints. Even the p-values relevant to the
  difference-in-slope parameters have been replaced by NA, since they are meaningless in
  this case, see \code{\link{davies.test}}.}
  \item{gap}{estimated coefficients, standard errors and t-values for the `gap' variables}
  \item{cov.var.diff}{if \code{var.diff=TRUE}, the covaraince matrix accounting for heteroscedastic errors.}
  \item{sigma.new}{if \code{var.diff=TRUE}, the square root of the estimated error variances in each interval.}
  \item{df.new}{if \code{var.diff=TRUE}, the residual degrees of freedom in each interval.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Vito M.R. Muggeo }
\seealso{ \code{\link{print.segmented}}, \code{\link{davies.test}} }
\examples{
##continues example from segmented()
# summary(segmented.model,short=TRUE)

## an heteroscedastic example..
# set.seed(123)
# n<-100
# x<-1:n/n
# y<- -x+1.5*pmax(x-.5,0)+rnorm(n,0,1)*ifelse(x<=.5,.4,.1)
# o<-lm(y~x)
# oseg<-segmented(o,seg.Z=~x,psi=.6)
# summary(oseg,var.diff=TRUE)$sigma.new
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
