\name{cub_spline_est}
\alias{cub_spline_est}
\title{
Cubic spline fitting 
}
\description{
The function cub_spline_est is an implementation of the (un)constrained cubic spline estimates proposed by Daouia, Noh and Park (2016).
}
\usage{
cub_spline_est(xtab, ytab, x, kn = ceiling((length(xtab))^(1/4)), method= "u",
               all.dea=FALSE, control = list("tm_limit" = 700))
}
\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{kn}{an integer specifying the number of inter-knot segments used in the computation of the spline estimate.}
  \item{method}{a character equal to "u" (unconstrained estimator), "m" (under the monotonicity constraint) or "mc" (under simultaneous monotonicity and concavity constraints).}
  \item{all.dea}{a boolean.}
  \item{control}{a list of parameters to the GLPK solver. See *Details* of help(Rglpk_solve_LP).}  
  }
\details{
Let \eqn{a} and \eqn{b} be, respectively, the minimum and maximum of the design points \eqn{x_1,\ldots,x_n}{x[1],...,x[n]}.  
Denote a partition of \eqn{[a,b]} by \eqn{a=t_0<t_1<\cdots<t_{k_n}=b}{a=t0<t1<...<tkn=b} (see below the selection process).
Let \eqn{N=k_n+3}{N=k[n]+3} and \eqn{\pi(x)=(\pi_0(x),\ldots,\pi_{N-1}(x))^T}{pi(x)=(pi[0](x),...,pi[N-1](x))^T} 
be the vector of normalized B-splines of order 4 based on the knot mesh \eqn{\{t_j\}}{t[j]} (see Daouia et al. (2015)). The
unconstrained (option \code{method="u"}) cubic spline estimate of the frontier \eqn{\varphi(x)}{varphi(x)} is then
defined by \eqn{\tilde\varphi_n(x)=\pi(x)^T \tilde\alpha}{varphihat[n](x)=pi(x)^T alphahat}, where \eqn{\tilde\alpha}{tilde(alpha)} minimizes
\deqn{\int_{0}^1\pi(x)^T \alpha \,dx = \sum_{j=1}^N \alpha_j \int_{0}^1\pi_j(x) \,dx}{integral from 0 to 1 pi(x)^T alpha dx = sum(from 1 to n alpha[j] integral from 0 to 1 pi(x)) dx}
over \eqn{\alpha\in\R^N}{alpha in R^n} subject to the envelopment constraints
\eqn{\pi(x_i)^T \alpha\geq y_i}{pi(x[i])^T alpha >= y[i]},  \eqn{i=1,\ldots,n}{i=1,...,n}.
A simple way of choosing the knot mesh in this unconstrained setting is by considering the  
\eqn{j/k_n}{j/k[n]}th quantiles \eqn{t_j = x_{[j n/k_n]}}{t[j] = x[[j*n/k[n]]]}  of the distinct values of \eqn{x_i}{x[i]} for \eqn{j=1,\ldots,k_n-1}{j=1,...,k[n]-1}. 
The number of inter-knot segments \eqn{k_n}{k[n]} is obtained by calling the function \code{\link{cub_spline_kn}} (option \code{method="u"}).

In what concerns the monotonicity constraint, we use the following suffcient condtion for the monotonicity:
\deqn{\alpha_0 \leq \alpha_1 \leq \cdots \leq \alpha_{N-1}}{alpha[0]<= alpha[1] ... <= alpha[N-1]}.
This condition was previously used in Lu et al. (2007) and Pya and Wood (2015). Note that since the condition corresponds to linear constraints on \eqn{\alpha}, the estimator satisfying the monotonocity can be obtained via linear programming.  

When the estimate is required to be both monotone and concave, we use the function \code{cub_spline_est} with the option \code{method="mc"}. Such estimate is obtained as the cubic spline function which minimizes the same linear objective function as the unconstrained estimate subject to the same linear envelopment constraints, the monotonicity constraint above and the additional linear concavity constraints \eqn{\pi''(t_j)^T \alpha\leq , j=0,1,\ldots,k_n}{pi''(t[j])^T alpha < 0, j=0,1,...,k[n]}, where the second derivative \eqn{\pi''}{pi''} is a linear spline. Regarding the choice of knots, we just apply the same scheme as for the unconstrained cubic spline estimate. 
 
}

\value{
Returns a numeric vector with the same length as \code{x}. Returns a vector of NA if no solution has been found by the solver (GLPK). 
}
\references{
Daouia, A., Noh, H. and Park, B.U. (2016). Data Envelope fitting with constrained polynomial splines. \emph{Journal of the Royal Statistical Society: Series B}, \bold{78}(1), 3-30. doi:10.1111/rssb.12098.

Lu, M., Zhang, Y. and Huang, J. (2007). Estimation of the mean function with panel count data using monotone polynomial splines. \emph{Biometrika}, \bold{94}, 705-718.

Pya, N. and Wood, S. N. (2015). Shape constrained additive models. \emph{Statistical Computing}, to appear.

Schumaker, L.L. (2007). \emph{Spline Functions: Basic Theory}, 3rd edition, Cambridge University Press.
}
\author{
Hohsuk Noh
}

\seealso{
\code{\link{cub_spline_kn}}
}
\examples{
data("air")
x.air <- seq(min(air$xtab), max(air$xtab), length.out=101)
 
# 1. Unconstrained cubic spline fits
# Optimal number of inter-knot segments via the BIC criterion
(kn.bic.air.u<-cub_spline_kn(air$xtab, air$ytab, 
 method="u", type="BIC"))
# Unconstrained spline estimate
y.cub.air.u<-cub_spline_est(air$xtab, air$ytab, 
 x.air, kn=kn.bic.air.u, method="u")

# 2. Monotonicity constraint
# Optimal number of inter-knot segments via the BIC criterion
(kn.bic.air.m<-cub_spline_kn(air$xtab, air$ytab, 
 method="m", type="BIC")) 
# Monotonic splines estimate
y.cub.air.m<-cub_spline_est(air$xtab, air$ytab, 
 x.air, kn=kn.bic.air.m, method="m")
   
# 3. Monotonicity and Concavity constraints
# Optimal number of inter-knot segments via the BIC criterion
(kn.bic.air.mc<-cub_spline_kn(air$xtab, air$ytab, 
 method="mc", type="BIC"))
# Monotonic/Concave splines estimate 
y.cub.air.mc<-cub_spline_est(air$xtab, air$ytab, 
 x.air, kn=kn.bic.air.mc, method="mc", all.dea=TRUE)

# Representation 
plot(x.air, y.cub.air.u, lty=1, lwd=4, col="green", 
 type="l", xlab="log(COST)", ylab="log(OUTPUT)")   
lines(x.air, y.cub.air.m, lty=2, lwd=4, col="cyan")
lines(x.air, y.cub.air.mc, lty=3, lwd=4, col="magenta")   
points(ytab~xtab, data=air)
legend("topleft", col=c("green","cyan","magenta"), 
lty=c(1,2,3), legend=c("unconstrained", "monotone", 
 "monotone + concave"), lwd=4, cex=0.8)    
}
\keyword{nonparametric}
\keyword{optimize}
