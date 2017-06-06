\name{poly_degree}
\alias{poly_degree}
\title{
AIC and BIC criteria for choosing the optimal degree of the polynomial frontier estimator
}
\description{
Computes the optimal degree of the unconstrained polynomial frontier estimator proposed by Hall, Park and Stern (1998).}
\usage{
poly_degree(xtab, ytab, prange=0:20, type="AIC", 
 control = list("tm_limit" = 700))
}
\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{prange}{a vector of integers specifying the range in which the optimal degree of the polynomial frontier estimator is to be selected.}
  \item{type}{a character equal to "AIC" or "BIC".}
  \item{control}{a list of parameters to the GLPK solver. See *Details* of help(Rglpk_solve_LP).}  
}
\details{
As the degree \eqn{p} of the polynomial estimator \eqn{\hat \varphi_{n,p}}{hat(varphi)[n,p]} (see \code{\link{poly_est}}) determines the dimensionality of the approximating function, we may view the problem of choosing p as model selection.
By analogy to the information criteria proposed by Daouia et al. (2015) in the boundary regression context, we obtain the optimal polynomial degree by minimizing
\deqn{
AIC(p) = \log \left( \sum_{i=1}^{n} (\hat \varphi_{n,p}(x_i)-y_i)\right) + (p+1)/n ,}{AIC(k) = log(sum from 1 to n (varphihat[n,p](x[i])-y[i])) + (p+1)/n,}
\deqn{BIC(p) = \log \left( \sum_{i=1}^{n} (\hat \varphi_{n,p}(x_i)-y_i)\right) + \log n (p+1)/(2n).}{AIC(k) = log(sum from 1 to n (varphihat[n,p](x[i])-y[i])) + log(n)(p+1)/(2n).}
The first one (option \code{type = "AIC"}) is similar to the famous Akaike information criterion Akaike (1973) and the second one
(option \code{type = "BIC"}) to the Bayesian information criterion Schwartz (1978). 
}
\value{
Returns an integer.
}
\references{
Akaike, H. (1973).  Information theory and an extension of the maximum likelihood principle, in \emph{Second International Symposium of Information Theory}, eds. B. N. Petrov and F. Csaki, Budapest: Akademia Kiado, 267--281.  

Daouia, A., Noh, H. and Park, B.U. (2015). Data Envelope fitting with constrained polynomial splines. \emph{Journal of the Royal Statistical Society: Series B}, to appear.

Hall, P., Park, B.U. and Stern, S.E. (1998). On polynomial estimators of frontiers and boundaries. \emph{Journal of Multivariate Analysis}, 66, 71-98.

Schwartz, G. (1978). Estimating the dimension of a model, \emph{Annals of Statistics}, 6, 461--464.
}
\author{
Hohsuk Noh.
}

\seealso{
\code{\link{poly_est}}
}
\examples{
data("air")
x.air <- seq(min(air$xtab), max(air$xtab), 
 length.out = 101)
# Optimal polynomial degrees via the AIC criterion
(p.aic.air <- poly_degree(air$xtab, air$ytab, 
 type = "AIC"))
# Optimal polynomial degrees via the BIC criterion  
(p.bic.air <- poly_degree(air$xtab, air$ytab, 
 type = "BIC"))
}
\keyword{nonparametric}
\keyword{optimize}
