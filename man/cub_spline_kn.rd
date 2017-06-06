\name{cub_spline_kn}
\alias{cub_spline_kn}
\title{
AIC and BIC criteria for choosing the number of inter-knot segments in cubic spline fits  
}
\description{
Computes the optimal number of inter-knot segments for the (un)constrained cubic spline fit proposed by Daouia, Noh and Park (2015).}
\usage{
cub_spline_kn(xtab, ytab, method, krange = 1:20, type = "AIC", 
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
The implementation of the unconstrained cubic spline smoother \eqn{\tilde\varphi_n}{tilde(varphi)[n]} (see \code{\link{cub_spline_est}}) 
is based on the knot mesh \eqn{\{t_j\}}{{t[j}}, with \eqn{t_j = x_{[j n/k_n]}}{t[j] = x[[j*n/k[n]]]} being the \eqn{j/k_n}{j/k[n]}th quantile 
of the distinct values of \eqn{x_i}{x[i]} for \eqn{j=1,\ldots,k_n-1}{j=1,...,k[n]-1}. 
Because the number of knots \eqn{k_n}{k[n]} determines the complexity of the spline approximation,
its choice may then be viewed as model selection through the minimization of the following two information criteria:
\deqn{
AIC(k) = \log \left( \sum_{i=1}^{n} (\tilde \varphi_n(x_i)- y_i) \right) + (k+2)/n,}{AIC(k) = log(sum from 1 to n (tilde(varphi)[n](x[i])-y[i])) + (k+2)/n,}
\deqn{BIC(k) = \log \left( \sum_{i=1}^{n} (\tilde \varphi_n(x_i) - y_i) \right) + \log n \cdot (k+2)/2n.}{BIC(k) = log(sum from 1 to n (tilde(varphi)[n](x[i])-y[i])) + log(n)(k+2)/2n.}
The first one (option \code{type = "AIC"}) is similar to the famous Akaike information criterion (Akaike, 1973) and the second one
(option \code{type = "BIC"}) to the Bayesian information criterion (Schwartz, 1978).
For the implementation of the concave cubic spline estimator, just apply the same scheme as for the unconstrained version.
}
\value{
Returns an integer.
}
\references{
Akaike, H. (1973).  Information theory and an extension of the maximum likelihood principle, in \emph{Second International Symposium of Information Theory}, eds. B. N. Petrov and F. Csaki, Budapest: Akademia Kiado, 267--281.  

Daouia, A., Noh, H. and Park, B.U. (2015). Data Envelope fitting with constrained polynomial splines. \emph{Journal of the Royal Statistical Society: Series B}, to appear.

Schwartz, G. (1978). Estimating the dimension of a model, \emph{Annals of Statistics}, 6, 461--464.
}
\author{
Hohsuk Noh.
}

\seealso{
\code{\link{cub_spline_est}}.
}
\examples{
data("air")
# a. Unconstrained cubic spline fits
(kn.bic.air.u<-cub_spline_kn(air$xtab, air$ytab, 
 method="u", type="BIC"))
# b. Monotone cubic spline smoother
(kn.bic.air.m<-cub_spline_kn(air$xtab, air$ytab, 
 method="m", type="BIC")) 
# c. Monotone and Concave cubic spline smoother
(kn.bic.air.mc<-cub_spline_kn(air$xtab, air$ytab, 
 method="mc", type="BIC"))
}
\keyword{nonparametric}
\keyword{optimize}
