\name{quad_spline_kn}
\alias{quad_spline_kn}
\title{
AIC and BIC criteria for choosing the optimal number of inter-knot segments in quadratic spline fits  
}
\description{
Computes the optimal number \eqn{k_n}{k[n]} of inter-knot segments in the quadratic spline fits proposed by Daouia, Noh and Park (2016).}
\usage{
quad_spline_kn(xtab, ytab, method, krange = 1:20, type = "AIC", 
 control = list("tm_limit" = 700))
}
\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{method}{a character equal to "u" (unconstrained estimator), "m" (under the monotonicity constraint) or "mc" (under simultaneous monotonicity and concavity constraints).}
  \item{krange}{a vector of integers specifying the range in which the optimal number of inter-knot segments is to be selected.}
  \item{type}{a character equal to "AIC" or "BIC".}
  \item{control}{a list of parameters to the GLPK solver. See *Details* of help(Rglpk_solve_LP).}  
}
\details{
For the implementation of the unconstrained quadratic spline smoother 
\eqn{\tilde\varphi_n}{tilde(varphi)[n]} (see \code{\link{quad_spline_est}}), based on the knot mesh 
\eqn{\{t_j = x_{[j n/k_n]}: j=1,\ldots,k_n-1\}}{{t[j] = x[[j*n/k[n]]], j=1,...,k[n]-1}},
the user has to employ the option \code{method="u"}. 
Since the number \eqn{k_n}{k[n]} determines the complexity of the spline approximation, 
its choice may be viewed as model selection via the minimization of the following Akaike (option \code{type="AIC"}) 
or Bayesian (option \code{type="BIC"}) information criteria:
\deqn{
A\tilde{I}C(k) = \log \left( \sum_{i=1}^{n} (\tilde \varphi_n(x_i)- y_i) \right) + (k+2)/n,}{AIC(k) = log(sum from 1 to n (tilde(varphi)[n](x[i])-y[i])) + (k+2)/n,}
\deqn{B\tilde{I}C(k) = \log \left( \sum_{i=1}^{n} (\tilde \varphi_n(x_i) - y_i) \right) + \log n \cdot (k+2)/2n.}{BIC(k) = log(sum from 1 to n (tilde(varphi)[n](x[i])-y[i])) + log(n)(k+2)/2n.}
For the implementation of the monotone (option \code{method="m"}) quadratic spline smoother \eqn{\hat\varphi_n}{hat(varphi)[n]} (see \code{\link{quad_spline_est}}),
the authors first suggest using the set of knots \eqn{\{ t_j = {\mathcal{X}_{[j \mathcal{N}/k_n]}},~j=1,\ldots,k_n-1 \}}{{t[j] = {X[jN/k[n]]}, j=1,...,k[n] }} 
among the FDH points \eqn{(\mathcal{X}_{\ell},\mathcal{Y}_{\ell})}{(X[l],Y[l])}, \eqn{\ell=1,\ldots,\mathcal{N}}{l=1,...,N}
(function \code{\link{quad_spline_est}}).
Then, they propose to choose \eqn{k_n}{k[n]} by minimizing the following AIC (option \code{type="AIC"}) or BIC (option \code{type="BIC"}) information criteria:
\deqn{
A\hat{I}C(k) = \log \left( \sum_{i=1}^{n} (\hat \varphi_n(x_i)- y_i) \right) + (k+2)/n,}{AIC(k) = log(sum from 1 to n (hat(varphi)[n](x[i])-y[i])) + (k+2)/n,}
\deqn{B\hat{I}C(k) = \log \left( \sum_{i=1}^{n} (\hat \varphi_n(x_i) - y_i) \right) + \log n \cdot (k+2)/2n.}{BIC(k) = log(sum from 1 to n (hat(varphi)[n](x[i])-y[i])) + log(n)(k+2)/2n.}
A small number of knots is typically needed as elucidated by the asymptotic theory. 

For the implementation of the monotone and concave (option \code{method="mc"}) spline estimator \eqn{\hat\varphi^{\star}_n}{varphihat*[n]}, 
just apply the same scheme as above by replacing the FDH points \eqn{(\mathcal{X}_{\ell},\mathcal{Y}_{\ell})}{(X[l],...,Y[l])}  
with the DEA points \eqn{(\mathcal{X}^*_{\ell},\mathcal{Y}^*_{\ell})}{(X*[l],...,Y*[l])}  (see \code{\link{dea_est}}). 
}
\value{
Returns an integer.
}
\references{
Akaike, H. (1973).  Information theory and an extension of the maximum likelihood principle, in \emph{Second International Symposium of Information Theory}, eds. B. N. Petrov and F. Csaki, Budapest: Akademia Kiado, 267--281.  

Daouia, A., Noh, H. and Park, B.U. (2016). Data Envelope fitting with constrained polynomial splines. \emph{Journal of the Royal Statistical Society: Series B}, \bold{78}(1), 3-30. doi:10.1111/rssb.12098.

Schwartz, G. (1978). Estimating the dimension of a model, \emph{Annals of Statistics}, 6, 461--464.
}
\author{
Hohsuk Noh.
}

\seealso{
\code{\link{quad_spline_est}}
}
\examples{
data("green")
\dontrun{
# BIC criteria for choosing the optimal number of 
# inter-knot segments in:   
# a. Unconstrained quadratic spline fits
(kn.bic.green.u <- quad_spline_kn(log(green$COST), 
 log(green$OUTPUT), method = "u", type = "BIC"))
# b. Monotone quadratic spline smoother
(kn.bic.green.m <- quad_spline_kn(log(green$COST), 
 log(green$OUTPUT), method = "m", type = "BIC"))  
# c. Monotone and concave quadratic spline smoother
(kn.bic.green.mc<-quad_spline_kn(log(green$COST), 
 log(green$OUTPUT), method = "mc", type = "BIC"))
}
}
\keyword{nonparametric}
\keyword{optimize}
