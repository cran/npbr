\name{poly_est}
\alias{poly_est}
\title{
Polynomial frontier estimators}
\description{
Computes the polynomial-type estimators of frontiers and boundaries proposed by Hall, Park and Stern (1998).
}
\usage{
poly_est(xtab, ytab, x, deg, control = list("tm_limit" = 700))
}
\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{deg}{an integer (polynomial degree).}
  \item{control}{a list of parameters to the GLPK solver. See *Details* of help(Rglpk_solve_LP).}  
}
\details{
The data edge is modeled by a single polynomial \eqn{\varphi_{\theta}(x) = \theta_0+\theta_1 x+\cdots+\theta_p x^p}{varphi[theta](x) = theta[0]+theta[1]x+...+theta[p]x^p} 
of known degree \eqn{p} that envelopes the full data and minimizes the area under its graph for \eqn{x\in[a,b]}{x in [a,b]}, with \eqn{a} and \eqn{b} 
being respectively the lower and upper endpoints of the design points \eqn{x_1,\ldots,x_n}{x[1],...,x[n]}.
The implemented function is the estimate \eqn{\hat\varphi_{n,p}(x) = \hat\theta_0+\hat\theta_1 x+\cdots+\hat\theta_p x^p}{varphihat[theta](x) = thetahat[0]+thetahat[1]x+...+thetahat[p]x^p}
 of \eqn{\varphi(x)}{varphi(x)}, where  \eqn{\hat\theta=(\hat\theta_0,\hat\theta_1,\cdots,\hat\theta_p)^T}{thetahat=(thetahat[0],thetahat[1],...,thetahat[p])^T}
 minimizes
\eqn{\int_{a}^b \varphi_{\theta}(x) \,dx}{integral from a to b phi[theta](x)dx}
over \eqn{\theta\in\R^{p+1}}{theta in R^(p+1)} subject to the envelopment constraints
\eqn{\varphi_{\theta}(x_i)\geq y_i}{varphi[theta](x[i]) >= y[i]},  \eqn{i=1,\ldots,n}{i=1,...,n}.
}
\value{
Returns a numeric vector with the same length as \code{x}. Returns a vector of NA if no solution has been found by the solver (GLPK). 
}


\references{
Hall, P., Park, B.U. and Stern, S.E. (1998). On polynomial estimators of frontiers and boundaries. \emph{Journal of Multivariate Analysis}, 66, 71-98. 
}
\author{
Hohsuk Noh.
}

\seealso{
\code{\link{loc_est}}
}
\examples{
data("air")
x.air <- seq(min(air$xtab), max(air$xtab), 
 length.out = 101)
# Optimal polynomial degrees via the AIC criterion
(p.aic.air <- poly_degree(air$xtab, air$ytab, 
 type = "AIC"))
# Polynomial boundaries estimate 
y.poly.air<-poly_est(air$xtab, air$ytab, x.air, 
 deg = p.aic.air)
# Representation
plot(x.air, y.poly.air, lty = 1, lwd = 4, 
 col = "magenta", type = "l")
points(ytab~xtab, data = air)  
legend("topleft",legend = paste("degree =", p.aic.air), 
 col = "magenta", lwd = 4, lty = 1)  
}
\keyword{optimize}