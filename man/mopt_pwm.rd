\name{mopt_pwm}
\alias{mopt_pwm}
\title{
Threshold selection for the PWM frontier estimator
}
\description{
This function implements the optimal smoothing parameter 
\code{coefm} involved in the probability-weighted moment frontier estimator of Daouia, Florens and Simar (2012).    
}
\usage{
mopt_pwm(xtab, ytab, x, a=2, rho, wind.coef=0.1)
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.} 
  \item{a}{a smoothing parameter (integer) larger than or equal to 2 (2 by default).}
  \item{rho}{a numeric vector of the same length as \code{x} or a scalar, which determines the values of rho.}  
  \item{wind.coef}{a scalar coefficient to be selected in the interval (0,1].}  
}

\details{
This is an implementation of an automated selection of the parameter \code{coefm}
involved in the probability-weighted moment (PWM) estimator \eqn{\tilde\varphi_{pwm}(x)}{tilde(varphi)[pwm](x)}  [see \code{\link{dfs_pwm}}].
It is an adaptation of the experimental method \code{\link{kopt_momt_pick}} by Daouia et al. (2010).
The idea is to select first (for each \eqn{x}) a grid of values for the parameter \code{coefm} given by 
\eqn{c = 1, \cdots, \min(10,[\sqrt{N_x}])}{c=1,...,min(10,sqrt(N[x]))}, where \eqn{N_x=\sum_{i=1}^n1_{\{x_i\le x\}}}{Nx=sum[i=1 to n]1(x[i]<=x)},  
and then select the \eqn{c} where the variation of the results is the smallest. 
To achieve this, we compute the standard deviations of \eqn{\tilde\varphi_{pwm}(x)}{tilde(varphi)[pwm](x)} over a ``window'' of size 
\eqn{wind.coef  \times \min(10,[\sqrt{N_x}])}{wind.coef x min(10,[sqrt(N[x])])}, where the coefficient \code{wind.coef} should be selected 
in the interval \eqn{(0,1]} in such a way to avoid numerical instabilities. 
The default option \code{wind.coef=0.1} corresponds to having a window large enough to cover around \eqn{10\%} 
of the possible values of \eqn{c} in the selected range of values for \code{coefm}. 
 The value of \eqn{c} where the standard deviation is minimal defines the desired \code{coefm}.   
}

\value{
Returns a numeric vector with the same length as \code{x}.
}

\references{
Daouia, A., Florens, J.-P. and Simar, L. (2010). Frontier estimation and extreme value theory. \emph{Bernoulli}, 16, 1039-1063.
}

\author{
Abdelaati Daouia and Thibault Laurent. 
}

\seealso{
\code{\link{dfs_pwm}}, \code{\link{kopt_momt_pick}}.
}


\examples{
data("post")
x.post<- seq(post$xinput[100],max(post$xinput), 
 length.out=100) 
\dontrun{
# When rho[x] is known and equal to 2:
best_cm.1<- mopt_pwm(post$xinput, post$yprod, 
 x.post, a=2, rho=2)
}
}

\keyword{nonparametric}