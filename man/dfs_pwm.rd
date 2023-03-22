\name{dfs_pwm}
\alias{dfs_pwm}
\title{
Probability-weighted moment frontier estimator
}
\description{
This function is an implementation of the probability-weighted moment frontier estimator developed by Daouia, Florens and Simar (2012).    
}
\usage{
dfs_pwm(xtab, ytab, x, coefm, a=2, rho, ci=TRUE)
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{coefm}{a tuning parameter (integer) larger than or equal to 1.}  
  \item{a}{a smoothing parameter (integer) larger than or equal to 2.}
  \item{rho}{a numeric vector of the same length as \code{x} or a scalar, which determines the values of rho.}
  \item{ci}{a boolean, TRUE for computing the confidence interval.} 
}
\details{
The regularized frontier estimator introduced by Daouia et al. (2012) is based 
on the unregularized probability-weighted moment estimator
\deqn{\hat{\varphi}_m(x) = \varphi_{fdh}(x) - \int_{0}^{\varphi_{fdh}(x)} \hat{F}^m(y|x)dy}{hat(varphi[m])(x)=varphi[fdh](x)-int[0 to varphi[fdh](x)]hat(F)^m(y|x)dy}
where the trimming order \eqn{m\geq 1}{m>=1} is an integer such that \eqn{m=m_n\to\infty}{m=m[n]->Infinity} as \eqn{n\to\infty}{n->Infinity}, 
and \eqn{\hat{F}(y|x)=\sum_{i=1}^n1_{(x_i\leq x,y_i\leq y)}/\sum_{i=1}^n1_{(x_i\leq x)}}{hat(F)(y|x)=sum[i=1]^n 1[(x[i] <= x, y[i] <=y)]/(sum[i=1]^n 1[(x[i] <= x])}.
The implemented estimator of \eqn{\varphi(x)}{varphi(x)} is then defined as 
\deqn{\tilde{\varphi}_m(x) = \hat{\varphi}_m(x) + \Gamma\left(1 + 1/\bar\rho_x\right)\left( 1/m\,\hat\ell_x\right)^{1/\bar\rho_x}}{tilde(varphi[m])(x) = hat(varphi[m])(x) + Gamma(1 + 1/hat(rho)[x])(1/(m.hat(l)[x]))^(1/hat(rho)[x])}
where
\deqn{\bar{\rho}_x    =   \log (a)\left\{ \log\Big( \frac{\hat\varphi_{m}(x)-\hat\varphi_{am}(x)}{\hat\varphi_{am}(x)-\hat\varphi_{a^2m}(x)}
\Big) \right\}^{-1} ,  \quad
\hat{\ell}_x     =  \frac {1}{m}\left[\frac{\Gamma(1+ 1/\bar\rho_x)\big(1-a^{-1/\bar\rho_x}\big)}{\hat\varphi_{m}(x)-\hat\varphi_{am}(x)}\right]^{\bar\rho_x},}{hat(rho)[x] = log(a){log((hat(varphi)[m](x)-hat(varphi)[am](x))/(hat(varphi)[am](x)-hat(varphi)[a^2m](x))}^{-1},
hat(l)[x] = 1/m ((Gamma(1+ 1/hat(rho)[x])(1-a^(-1/hat(rho)[x])))/(hat(varphi)[m](x)-hat(varphi)[am](x)))^(hat(rho)[x]),}
with \eqn{a\geq 2}{a >= 2} being a fixed integer. If the true tail-index \eqn{\rho_x=\beta_x+2}{rho[x]=beta[x]+2} is known, we set \eqn{\bar{\rho}_x=\rho_x}{hat(rho)[x]=rho[x]} in the expressions above.
The two smoothing parameters \eqn{m} and \eqn{a} have to be fixed by the user in the 4th and 5th arguments of the function.

The pointwise \eqn{95\%} confidence interval of \eqn{\varphi(x)}{varphi(x)} derived from the asymptotic normality of \eqn{\tilde\varphi_{m}(x)}{tilde(varphi)[m](x)} 
is given by \eqn{[\tilde{\varphi}_m(x)  \pm   1.96 \, \hat\sigma(m,x)/\sqrt{n}]}{[tilde(varphi)[m](x) +- 1.96 hat(sigma)(m,x)/sqrt(n)]}
where 
\deqn{ 
 \hat\sigma^2(m,x)= \frac{2m^2}{ \hat F_X(x)}\int_0^{\varphi_{fdh}(x)}  
\int_0^{\varphi_{fdh}(x)}  \hat F^{m}(y|x)\hat F^{m-1}(u|x)(1-\hat F(u|x))
1_{(y\le u)}\, dy\,du,}{hat(sigma)^2(m,x)= (2m^2)/( hat(F)[X](x))int[0 -> varphi[fdh](x)]int[0 -> varphi[fdh](x)] hat(F)^m(y|x)hat(F)^(m-1)(u|x)(1-hat(F)(u|x)) 1[(y <= u)]dydu,}
with \eqn{\hat F_X(x) =(1/n)\sum_{i=1}^n1_{(x_i\leq x)}}{hat(F)[X](x) =(1/n)sum[i=1 -> n] 1[(x[i] <= x)]}. 
Note that the standard deviation \eqn{\sigma(m,x)/\sqrt{n}}{sigma(m,x)/sqrt(n)} of the bias-corrected estimator \eqn{\tilde{\varphi}_m(x)}{tilde(varphi)[m](x)} is adjusted by a bootstrap estimator 
in the numerical illustrations of Daouia et al. (2012), whereas the exact estimate \eqn{\hat\sigma(m,x)/\sqrt{n}}{hat(sigma)(m,x)/sqrt(n)} is utilized in the implemented function.
A practical choice of \eqn{m} that Daouia et al. (2012) have employed is the simple rule of thumb 
\eqn{m=coefm \times N^{1/3}_x}{m=coefm N^(1/3)[x]}, where \eqn{N_x=\sum_{i=1}^n1_{\{x_i\le x\}}}{N[x]=sum[i=1 -> n] 1[{x[i] <= x}]}, 
and the integer \code{coefm} as well as the second smoothing parameter \code{a} are to be tuned by the user to avoid numerical instabilities 
in the pointwise estimates of the tail-index \eqn{\rho_x}{rho[x]} and the frontier function \eqn{\varphi(x)}{varphi(x)}. 
The user may start with the values \code{coefm=5} and \code{a=2} [respectively, \code{coefm=10} and \code{a=20}] 
for computing the estimator \eqn{\tilde{\varphi}_m(x)}{tilde(varphi)[m](x)}  [respectively, \eqn{\bar{\rho}_x}{hat(rho)[x]}]. Note that tail-index estimation and frontier estimation are conducted separately. 
}

\value{
Returns a numeric vector with the same length as \code{x}.
}


\note{
The computational burden here is demanding, so be forewarned. 
}

\references{
Daouia, A., Florens, J.-P. and Simar, L. (2012). Regularization of Nonparametric Frontier Estimators. \emph{Journal of Econometrics}, 168, 285-299.
}

\author{
Abdelaati Daouia and Thibault Laurent (converted from Abdelaati Daouia's Matlab code). 
}

\seealso{
\code{\link{rho_pwm}}, \code{\link{mopt_pwm}}.
}


\examples{
data("post")
x.post<- seq(post$xinput[100],max(post$xinput), 
 length.out=100) 
\dontrun{
# 1. When rho[x] is known and equal to 2, we set:
rho<-2
res.pwm.1<- dfs_pwm(post$xinput, post$yprod, x.post, coefm=5,
 a=2, rho, ci=TRUE)
# 2. When rho[x] is unknown and dependent of x, 
# its estimate hat(rho[x]) is obtained via:
rho_pwm <- rho_pwm(post$xinput, post$yprod, x.post, coefm=10, a=20)
# and the corresponding frontier estimator via: 
res.pwm.2<- dfs_pwm(post$xinput, post$yprod, x.post, coefm=5,
 a=2, rho_pwm, ci=TRUE)
# 3. When rho[x] is unknown but independent of x, 
# a robust estimation strategy is by using the (trimmed) mean 
# over the estimates hat(rho[x]): 
rho_trimmean<-mean(rho_pwm, trim=0.00)
res.pwm.3<- dfs_pwm(post$xinput, post$yprod, x.post, coefm=5,
 a=2, rho_trimmean, ci=TRUE)
}
}

\keyword{nonparametric}