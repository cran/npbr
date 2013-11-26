\name{loc_est}
\alias{loc_est}
\title{
Local linear frontier estimator
}
\description{
Computes the local linear smoothing frontier estimator of Hall, Park and Stern (1998).    
}
\usage{
loc_est(xtab, ytab, x, h)
}

\arguments{
  \item{xtab}{numeric vectors containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}}
  \item{ytab}{numeric vectors of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed}
  \item{h}{determines the bandwidth at which the local linear estimate will be computed}
}
\details{
The implemented local linear estimator of \eqn{\varphi(x)}{varphi(x)} is defined by
\deqn{
\hat \varphi_{n,LL}(x) = \min  \Big\{ z : {\rm there~exists~} \theta\geq 0 ~{\rm such~that~} y_i \leq z + \theta(x_i - x)}{hat(varphi)[n,LL](x) = min{ z : there exists  theta>=0 such that y_i <= z + theta(xi - x)}}
\deqn{{\rm for~all}~i~{\rm such~that~}x_i \in (x-h,x+h) \Big\}.}{for all i such that xi in (x-h,x+h).}
Hall and Park (2004) proposed a bootstrap procedure for selecting the optimal bandwidth \eqn{h} in  \eqn{\hat \varphi_{n,LL}}{hatvarphi[n,LL]}. See the function \code{\link{loc_est_bw}}.
}
\value{
Returns a numeric vector with the same length as \code{x}
}
\references{
Hall, P. and Park, B.U. (2004). Bandwidth choice for local polynomial estimation of smooth boundaries. \emph{Journal of Multivariate Analysis}, 91 (2), 240-261. 

Hall, P., Park, B.U. and Stern, S.E. (1998). On polynomial estimators of frontiers and boundaries. \emph{Journal of Multivariate Analysis}, 66, 71-98.  
}

\author{
Hohsuk Noh
}

\seealso{
\code{\link{poly_est}}, \code{\link{loc_est_bw}}
}


\examples{
data("nuclear")
plot(ytab~xtab, data=nuclear)
x <- seq(min(nuclear$xtab), max(nuclear$xtab), length.out=1001)
lines(x, loc_est(nuclear$xtab, nuclear$ytab, x, h=79.12), lty=1, col="red")   
lines(x, loc_est(nuclear$xtab, nuclear$ytab, x, h=40), lty=2, col="blue") 
legend("topleft",legend=c("h=72.70", "h=40"), col=c("red","blue"), 
lty=c(1,2))
}

\keyword{optimize}