\name{rho_pwm}
\alias{rho_pwm}
\title{
Probability-weighted moment frontier estimator
}
\description{
This function is an implementation of the Probability-weighted moment frontier estimator developed by Daouia, Florens and Simar (2012).    
}
\usage{
rho_pwm(xtab, ytab, x, a=2, lrho=1, urho=Inf)
}
\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{a}{a smoothing parameter (integer) larger than or equal to 2.}
  \item{lrho}{a scalar, minimum rho threshold value.}  
  \item{urho}{a scalar, maximum rho threshold value.}  
}

\details{
The function computes the probability-weighted moment (PWM) estimator \eqn{\bar\rho_x}{bar(rho)[x]} utilized in the frontier estimate 
\eqn{\tilde\varphi_{pwm}(x)}{tilde(varphi)[pwm](x)}[see \code{\link{dfs_pwm}}]. 
This estimator depends on the smoothing parameters \eqn{a} and \eqn{m}. A simple selection rule of thumb that Daouia et al. (2012) have employed is 
\eqn{a=2} 
[default option in the 4th argument of the function] 
and \eqn{m=coefm \times N^{1/3}_x}{m=coefm x N[x]^(1/3)}, where \eqn{N_x=\sum_{i=1}^n1_{\{x_i\le x\}}}{Nx=sum[i=1 to n]1(x[i]<=x)} 
and the integer \code{coefm} is to be tuned by the user. 
To choose this parameter in an optimal way for each \eqn{x}, we adapt the automated threshold selection method of Daouia et al. (2010) as follows:  
We first evaluate the estimator \eqn{\bar\rho_x}{bar(rho)[x]} over a grid of values of \code{coefm} given by 
\eqn{c = 1, \cdots, 150}{c=1,...,150}. 
Then, we select the \eqn{c} where the variation of the results is the smallest. This is achieved by computing the standard deviation of  the estimates \eqn{\bar\rho_x}{bar(rho)[x]} over a ``window'' of 
\eqn{\max([\sqrt{150}],3)}{max(sqrt(150),3)} successive values of \eqn{c}. The value of \eqn{c} where this standard deviation is minimal defines the value of \code{coefm}.
The user can also appreciably improve the estimation of the extreme-value index \eqn{\rho_x}{rho[x]} and the frontier function \eqn{\varphi_x}{varphi[x]} itself by tuning the choice of the lower limit 
(default option \code{lrho=1}) and upper limit (default option \code{urho=Inf}).
}

\value{
Returns a numeric vector with the same length as \code{x}.
}

\note{
The computational burden here is demanding, so be forewarned. 
}

\references{
Daouia, A., Florens, J.-P. and Simar, L. (2010). Frontier estimation and extreme value theory. \emph{Bernoulli}, 16, 1039-1063.

Daouia, A., Florens, J.-P. and Simar, L. (2012). Regularization of Nonparametric Frontier Estimators. \emph{Journal of Econometrics}, 168, 285-299.
}

\author{
Abdelaati Daouia and Thibault Laurent. 
}

\seealso{
\code{\link{dfs_pwm}}, \code{\link{mopt_pwm}}.
}


\examples{
data("post")
x.post<- seq(post$xinput[100],max(post$xinput), 
 length.out=100) 
\dontrun{
# When rho[x] is unknown and dependent of x, 
# its estimate hat(rho[x]) is obtained via:
rho_pwm <- rho_pwm(post$xinput, post$yprod, x.post,  a=20)
}
}

\keyword{nonparametric}