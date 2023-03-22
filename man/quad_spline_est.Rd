\name{quad_spline_est}
\alias{quad_spline_est}
\title{
Quadratic spline frontiers 
}
\description{
This function is an implementation of the (un)constrained quadratic spline smoother proposed by Daouia, Noh and Park (2016).
}
\usage{
quad_spline_est(xtab, ytab, x, kn = ceiling((length(xtab))^(1/4)), method= "u", 
 all.dea = FALSE, control = list("tm_limit" = 700))
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
Let \eqn{N=k_n+1}{N=k[n]+2} and \eqn{\pi(x)=(\pi_1(x),\ldots,\pi_N(x))^T}{pi(x)=(pi[1](x),...,pi[N](x))^T} be the vector of normalized 
B-splines of order 3 based on the knot mesh \eqn{\{t_j\}}{t[j]} (see, \emph{e.g.}, Schumaker (2007)).
When the true frontier \eqn{\varphi(x)}{varphi(x)} is known or required to be monotone nondecreasing (option \code{cv=0}), 
its constrained quadratic spline estimate is defined by  \eqn{\hat\varphi_n(x)=\pi(x)^T \hat\alpha}{varphihat[n](x)=pi(x)^T alphahat}, where \eqn{\hat\alpha}{alphahat} minimizes
\deqn{\int_{0}^1\pi(x)^T \alpha \,dx = \sum_{j=1}^N \alpha_j \int_{0}^1\pi_j(x) \,dx}{integral from 0 to 1 pi(x)^T alpha dx = sum(from 1 to n alpha[j] integral from 0 to 1 pi(x)) dx}
over \eqn{\alpha\in\R^N}{alpha in R^n} subject to the envelopment and monotonicity constraints
\eqn{\pi(x_i)^T \alpha\geq y_i}{pi(x[i])^T alpha >= y[i]},  \eqn{i=1,\ldots,n}{i=1,...,n}, and \eqn{\pi'(t_j)^T \alpha\geq 0}{pi'(t[j])^T alpha >= 0}, \eqn{j=0,1,\ldots,k_n}{j=0,1,...,k[n]}, 
with \eqn{\pi'}{pi'} being the derivative of \eqn{\pi}{pi}.

Considering the special connection of the spline smoother \eqn{\hat \varphi_n}{varphihat[n]} with the traditional FDH frontier
\eqn{\varphi_n}{varphi[n]} (see the function \code{\link{dea_est}}),
Daouia et al. (2015) propose an easy way of choosing the knot mesh.
Let \eqn{(\mathcal{X}_1,\mathcal{Y}_1),\ldots, (\mathcal{X}_\mathcal{N},\mathcal{Y}_\mathcal{N})}{(X[1],Y[1]),...,(X[N],Y[N])} 
be the observations \eqn{(x_i,y_i)}{(x[i],y[i])} lying on the FDH boundary (\emph{i.e.} \eqn{y_i=\varphi_n(x_i)}{y[i]=varphi[n](x[i])}).
The basic idea is to pick out a set of knots equally spaced in percentile ranks 
among the \eqn{\mathcal{N}}{N} FDH points \eqn{(\mathcal{X}_{\ell},\mathcal{Y}_{\ell})}{(X[l],Y[l])} by
taking \eqn{t_j = {\mathcal{X}_{[j \mathcal{N}/k_n]}}}{t[j] = {X}[jN/k[n]]}, the \eqn{j/k_n}{j/k[n]}th quantile 
of the values of \eqn{\mathcal{X}_{\ell}}{X[l]} for \eqn{j=1,\ldots,k_n-1}{j=1,...,k[n]-1}.
The choice of the number of internal knots is then viewed as model selection 
through the minimization of the AIC and BIC information criteria (see the function \code{\link{quad_spline_kn}}).

When the monotone boundary \eqn{\varphi(x)}{varphi(x)}  is also believed to be concave (option \code{cv=1}),
its constrained fit is defined as \eqn{\hat\varphi^{\star}_n(x)=\pi(x)^T \hat\alpha^{\star}}{varphihat*[n](x)=pi(x)^T alphahat*}, where \eqn{\hat\alpha^{\star}\in\R^N}{alphahat* in R^N} 
minimizes the same objective function as \eqn{\hat\alpha}{alphahat} subject to the same envelopment 
and monotonicity constraints and the additional concavity constraints
\eqn{\pi''(t^*_j)^T \alpha\leq 0}{pi''(t^*[j]^T alpha <= 0}, \eqn{j=1,\ldots,k_n,}{j=1,...,k[n]}
where \eqn{\pi''}{pi''} is the constant second derivative of \eqn{\pi}{pi} on each inter-knot interval  and \eqn{t^*_j}{t*[j]} is the midpoint of \eqn{(t_{j-1},t_j]}{(t[j-1],t[j])}.

Regarding the choice of knots, the same scheme as for \eqn{\hat\varphi_n}{varphihat[n]} can be applied by replacing 
the FDH points \eqn{(\mathcal{X}_1,\mathcal{Y}_1),\ldots, (\mathcal{X}_\mathcal{N},\mathcal{Y}_\mathcal{N})}{(X[1],Y[1]),...,(X[N],Y[N])}   
with the DEA points \eqn{(\mathcal{X}^*_1,\mathcal{Y}^*_1),\ldots, (\mathcal{X}^*_\mathcal{M},\mathcal{Y}^*_\mathcal{M})}{(X*[1],Y*[1]),...,(X*[M],Y*[M])}, that is, 
the observations \eqn{(x_i,y_i)}{(x[i],y[i])} lying on the piecewise linear DEA frontier  (see the function \code{\link{dea_est}}). 
Alternatively, the strategy of just using all the DEA points as knots is also 
working quite well for datasets of modest size as shown in Daouia et al. (2016). 
In this case, the user has to choose the option \code{all.dea=TRUE}.
}

\value{
Returns a numeric vector with the same length as \code{x}. Returns a vector of NA if no solution has been found by the solver (GLPK). 
}

\references{
Daouia, A., Noh, H. and Park, B.U. (2016). Data Envelope fitting with constrained polynomial splines. \emph{Journal of the Royal Statistical Society: Series B}, \bold{78}(1), 3-30. doi:10.1111/rssb.12098.

Schumaker, L.L. (2007). \emph{Spline Functions: Basic Theory}, 3rd edition, Cambridge University Press.
}
\author{
Hohsuk Noh.
}

\seealso{
\code{\link{quad_spline_kn}}
}
\examples{
\dontrun{
data("green")
x.green <- seq(min(log(green$COST)), max(log(green$COST)), length.out=101)
# 1. Unconstrained quadratic spline fits
# Optimal number of inter-knot segments via the BIC criterion
(kn.bic.green.u<-quad_spline_kn(log(green$COST), 
 log(green$OUTPUT), method="u", type="BIC"))
# Unconstrained spline estimate
y.quad.green.u<-quad_spline_est(log(green$COST), 
 log(green$OUTPUT), x.green, kn=kn.bic.green.u, method="u")

# 2. Monotonicity constraint
# Optimal number of inter-knot segments via the BIC criterion
(kn.bic.green.m<-quad_spline_kn(log(green$COST), 
 log(green$OUTPUT), method="m", type="BIC"))
# Monotonic splines estimate
y.quad.green.m<-quad_spline_est(log(green$COST), 
 log(green$OUTPUT), x.green, kn=kn.bic.green.m, method="m") 
   
# 3. Monotonicity and Concavity constraints
# Optimal number of inter-knot segments via the BIC criterion
(kn.bic.green.mc<-quad_spline_kn(log(green$COST), 
 log(green$OUTPUT), method="mc", type="BIC"))
# Monotonic/Concave splines estimate 
y.quad.green.mc<-quad_spline_est(log(green$COST), 
 log(green$OUTPUT), x.green, kn=kn.bic.green.mc, 
 method="mc", all.dea=TRUE)

# Representation 
plot(x.green, y.quad.green.u, lty=1, lwd=4, col="green", 
 type="l", xlab="log(COST)", ylab="log(OUTPUT)")   
lines(x.green, y.quad.green.m, lty=2, lwd=4, col="cyan")
lines(x.green, y.quad.green.mc, lwd=4, lty=3, col="magenta")   
points(log(OUTPUT)~log(COST), data=green)
legend("topleft", col=c("green","cyan","magenta"), 
lty=c(1,2,3), legend=c("unconstrained", "monotone", 
 "monotone + concave"), lwd=4, cex=0.8) 
}   
}
\keyword{nonparametric}
\keyword{optimize}
