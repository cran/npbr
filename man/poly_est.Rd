\name{poly_est}
\alias{poly_est}
\title{
Polynomial-based edge estimators}
\description{
Computes polynomial-type estimators of frontiers and boundaries (Hall et al., 1998)
}
\usage{
poly_est(xtab, ytab, x, deg)
}
\arguments{
  \item{xtab}{numeric vectors containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}}
  \item{ytab}{numeric vectors of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed}
  \item{deg}{an integer (polynomial degree)}
}
\details{
The data edge is modeled by a single polynomial \eqn{\varphi_{\theta}(x) = \theta_0+\theta_1 x+\cdots+\theta_p x^p}{varphi[theta](x) = theta[0]+theta[1]x+...+theta[p]x^p} 
of known degree \eqn{p} that envelopes the full data and minimizes the area under its graph for \eqn{x\in[a,b]}{x in [a,b]}, with \eqn{a} and \eqn{b} 
being respectively the lower and upper endpoints of the design points \eqn{x_1,\ldots,x_n}{x[1],...,x[n]}.
The function is the estimate \eqn{\hat\varphi_{n,P}(x) = \hat\theta_0+\hat\theta_1 x+\cdots+\hat\theta_p x^p}{varphihat[theta](x) = thetahat[0]+thetahat[1]x+...+thetahat[p]x^p}
 of \eqn{\varphi(x)}{varphi(x)}, where  \eqn{\hat\theta=(\hat\theta_0,\hat\theta_1,\cdots,\hat\theta_p)^T}{thetahat=(thetahat[0],thetahat[1],...,thetahat[p])^T}
 minimizes
\eqn{\int_{a}^b \varphi_{\theta}(x) \,dx}{integral from a to b phi[theta](x)dx}
over \eqn{\theta\in\R^{p+1}}{theta in R^(p+1)} subject to the envelopment constraints
\eqn{\varphi_{\theta}(x_i)\geq y_i}{varphi[theta](x[i]) >= y[i]},  \eqn{i=1,\ldots,n}{i=1,...,n}.
}
\value{
Returns a vector of numeric with the same length than \code{x}
}

\references{
Hall, P., Park, B.U. and Stern, S.E. (1998). On polynomial estimators of frontiers and boundaries. \emph{Journal of Multivariate Analysis}, 66, 71-98. 
}
\author{
Hohsuk Noh
}

\seealso{
\code{\link{loc_est}}
}
\examples{
data("nuclear")
plot(ytab~xtab, data=nuclear)
x <- seq(min(nuclear$xtab), max(nuclear$xtab), length.out=1001)
lines(x, poly_est(nuclear$xtab, nuclear$ytab, x, deg=2), lty=1, col="red")   
lines(x, poly_est(nuclear$xtab, nuclear$ytab, x, deg=4), lty=2, col="blue") 
legend("topleft",legend=c("degree=2", "degree=4"), col=c("red","blue"), 
lty=c(1,2))
}
\keyword{optimize}