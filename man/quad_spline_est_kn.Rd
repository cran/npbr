\name{quad_spline_est_kn}
\alias{quad_spline_est_kn}
\title{
AIC (or BIC) criterion for choosing the number of knots for quadratic splines  
}
\description{
Computes the optimal number of knots for the constrained quadratic spline fit proposed by Daouia, Noh and Park (2013)}
\usage{
quad_spline_est_kn(xtab, ytab, x, cv, krange = 1:20, type = "AIC")
}
\arguments{
  \item{xtab}{numeric vectors containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}}
  \item{ytab}{numeric vectors of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed}
  \item{cv}{an integer equal to 0 (constraint of monotonicity only) or 1 (both constrainst of monotonicity and concavity)}
  \item{krange}{a vector of integer specifying the number of knots at which the spline estimate will be computed}
  \item{type}{a character equal to "AIC" or "BIC"}
}
\details{
For the implementation of the monotone quadratic spline smoother \eqn{\hat\varphi_n}{varphihat[n]},
Daouia \emph{et al.} (2013) first suggest using the set of knots \eqn{\{ t_j = {\mathcal{X}_{[j \mathcal{N}/k_n]}},~j=1,\ldots,k_n-1 \}}{{t[j] = {X[jN/k[n]]}, j=1,...,k[n] }} 
among the FDH points \eqn{(\mathcal{X}_{\ell},\mathcal{Y}_{\ell})}{(X[l],Y[l])}, \eqn{\ell=1,\ldots,\mathcal{N}}{l=1,...,N}
(function \code{\link{quad_spline_est}}).
Because the number of knots \eqn{k_n}{k[n]} determines the complexity of the spline approximation,
its choice may then be viewed as model selection
 through the minimization of the following two information criteria:
\deqn{
AIC(k) = \log \left( \sum_{i=1}^{n} |y_i- \hat \varphi_n(x_i)|\right) + 2(k+2)/n,}{AIC(k) = log(sum from 1 to n |y[i]- varphihat[n](x[i])|) + 2(k+2)/n,}
\deqn{BIC(k) = \log \left( \sum_{i=1}^{n} |y_i- \hat \varphi_n(x_i)|\right) + \log n \cdot (k+2)/n.}{BIC(k) = log(sum from 1 to n |y[i]- varphihat[n](x[i])|) + log(n)(k+2)/n.}
The first one (option \code{type = "AIC"}) is similar to the famous Akaike information criterion (Akaike, 1973) and the second one
(option \code{type = "BIC"}) to the Bayesian information criterion (Schwartz, 1978).
A small number of knots is typically needed as elucidated by the asymptotic theory.

For the implementation of the monotone and concave spline estimator \eqn{\hat\varphi^{\star}_n}{varphihat*[n]}, 
just apply the same scheme as above by replacing the FDH points \eqn{(\mathcal{X}_{\ell},\mathcal{Y}_{\ell})}{(X[l],...,Y[l])} 
with the DEA points \eqn{(\mathcal{X}^*_{\ell},\mathcal{Y}^*_{\ell})}{(X*[l],...,Y*[l])} (see \code{\link{dea_est}}). 
}
\value{
Returns an integer
}
\references{
Akaike, H. (1973).  Information theory and an extension of the maximum likelihood principle, in \emph{Second International Symposium of Information Theory}, eds. B. N. Petrov and F. Csaki, Budapest: Akademia Kiado, 267--281.  

Daouia, A., Noh, H. and Park, B.U. (2013). Data Envelope Fitting with Constrained Polynomial splines. \emph{TSE Working Paper}, \url{http://www.tse-fr.eu/images/doc/wp/etrie/wp_tse_449.pdf}.

Schwartz, G. (1978). Estimating the dimension of a model, \emph{Annals of Statistics}, 6, 461--464.
}
\author{
Hohsuk Noh
}

\seealso{
\code{\link{quad_spline_est}}
}
\examples{
data("green")
x <- seq(min(log(green$COST)), max(log(green$COST)), length.out=1001)
quad_spline_est_kn(log(green$COST), log(green$OUTPUT), x, cv=1, type="AIC")   
quad_spline_est_kn(log(green$COST), log(green$OUTPUT), x, cv=1, type="BIC")
}
\keyword{nonparametric}
\keyword{optimize}
