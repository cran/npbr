\name{loc_est}
\alias{loc_est}
\title{
Local linear frontier estimator
}
\description{
Computes the local linear smoothing frontier estimator of Hall, Park and Stern (1998) and Hall and Park (2004).    
}
\usage{
loc_est(xtab, ytab, x, h, method="u", control = list("tm_limit" = 700))
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{h}{determines the bandwidth at which the local linear estimate will be computed.}
  \item{method}{a character equal to "u" (unconstrained estimator) or "m" (improved version of the unconstrained estimator).} 
  \item{control}{a list of parameters to the GLPK solver. See *Details* of help(Rglpk_solve_LP).}  
}
\details{
In the unconstrained case (option \code{method="u"}), the implemented estimator of \eqn{\varphi(x)}{varphi(x)} is defined by
\deqn{
\hat \varphi_{n,LL}(x) = \min  \Big\{ z : {\rm there~exists~} \theta ~{\rm such~that~} y_i \leq z + \theta(x_i - x)}{hat(varphi)[n,LL](x) = min{ z : there exists  theta such that y_i <= z + theta(xi - x)}}
\deqn{{\rm for~all}~i~{\rm such~that~}x_i \in (x-h,x+h) \Big\},}{for all i such that xi in (x-h,x+h),}
where the bandwidth \eqn{h} has to be fixed by the user in the 4th argument of the function.
This estimator may lack of smoothness in case of small samples and has no guarantee of being monotone  even if the true frontier is so.  
Following the curvature of the monotone frontier \eqn{\varphi}{varphi}, the unconstrained  estimator \eqn{\hat \varphi_{n,LL}}{hat(varphi)[n,LL]}  is likely to exhibit substantial bias, especially at the sample boundaries 
(see Daouia et al (2016) for numerical illustrations). A simple way to remedy to this drawback is by imposing the extra condition \eqn{\theta \geq 0}{theta <0} in the definition of \eqn{\hat \varphi_{n,LL}(x)}{hat(varphi)[n,LL](x)} to get 
\deqn{
\tilde \varphi_{n,LL}(x) = \min  \Big\{ z : {\rm there~exists~} \theta\geq 0 ~{\rm such~that~} y_i \leq z + \theta(x_i - x)}{hat(varphi)[n,LL](x) = min{ z : there exists  theta>=0 such that y_i <= z + theta(xi - x)}}
\deqn{{\rm for~all}~i~{\rm such~that~}x_i \in (x-h,x+h) \Big\}.}{for all i such that xi in (x-h,x+h).}
As shown in Daouia et al (2016), this version only reduces the vexing bias and border defects of the original estimator when the true frontier is monotone. 
The option \code{method="m"} indicates that the improved fit \eqn{\tilde \varphi_{n,LL}(x)}{tilde(varphi)[n,LL](x)} should be utilized in place of \eqn{\hat \varphi_{n,LL}(x)}{hat(varphi)[n,LL](x)}. 
Hall and Park (2004) proposed a bootstrap procedure for selecting the optimal 
bandwidth \eqn{h} in  \eqn{\hat \varphi_{n,LL}(x)}{hat(varphi)[n,LL](x)} and \eqn{\tilde \varphi_{n,LL}(x)}{tilde(varphi)[n,LL](x)} (see the function \code{\link{loc_est_bw}}).
}

\value{
Returns a numeric vector with the same length as \code{x}. Returns a vector of NA if no solution has been found by the solver (GLPK). 
}

\references{
Daouia, A., Noh, H. and Park, B.U. (2016). Data Envelope fitting with constrained polynomial splines. \emph{Journal of the Royal Statistical Society: Series B}, \bold{78}(1), 3-30. doi:10.1111/rssb.12098.

Hall, P. and Park, B.U. (2004). Bandwidth choice for local polynomial estimation of smooth boundaries. \emph{Journal of Multivariate Analysis}, 91, 240-261. 

Hall, P., Park, B.U. and Stern, S.E. (1998). On polynomial estimators of frontiers and boundaries. \emph{Journal of Multivariate Analysis}, 66, 71-98.  
}

\author{
Hohsuk Noh.
}

\seealso{
\code{\link{loc_est_bw}}.
}


\examples{
data("nuclear")
x.nucl <- seq(min(nuclear$xtab), max(nuclear$xtab), 
 length.out=101) 
# 1. Unconstrained estimator
# Optimal bandwidths over 100 bootstrap replications
\dontrun{
h.nucl.u <- loc_est_bw(nuclear$xtab, nuclear$ytab, 
 x.nucl, h=40, B=100, method="u")
}
(h.nucl.u<-79.11877)
y.nucl.u<-loc_est(nuclear$xtab, nuclear$ytab, x.nucl, 
 h=h.nucl.u, method="u")

# 2. improved version of the estimator
# Optimal bandwidths over 100 bootstrap replications
\dontrun{ 
h.nucl.m <- loc_est_bw(nuclear$xtab, nuclear$ytab, 
 x.nucl, h=40, B=100, method="m")
} 
(h.nucl.m<-79.12)
y.nucl.m<-loc_est(nuclear$xtab, nuclear$ytab, x.nucl, 
 h=h.nucl.m, method="m") 

# 3. Representation 
plot(x.nucl, y.nucl.u, lty=1, lwd=4, col="magenta", type="l")
lines(x.nucl, y.nucl.m, lty=2, lwd=4, col="cyan") 
points(ytab~xtab, data=nuclear)
legend("topleft",legend=c("unconstrained", "improved"), 
 col=c("magenta","cyan"), lwd=4, lty=c(1,2))
}

\keyword{optimize}
