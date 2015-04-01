\name{dfs_momt}
\alias{dfs_momt}
\title{
Moment frontier estimator
}
\description{
This function is an implementation of the moment-type estimator developed by Daouia, Florens and Simar (2010).    
}
\usage{
dfs_momt(xtab, ytab, x, rho, k, ci=TRUE)
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{rho}{a numeric vector of the same length as \code{x} or a scalar, which determines the values of rho.}  
  \item{k}{a numeric vector of the same length as \code{x} or a scalar, which determines the thresholds at which the moment estimator will be computed.}
  \item{ci}{a boolean, TRUE for computing the confidence interval.} 
}
\details{
Combining ideas from Dekkers, Einmahl and de Haan (1989) with the dimensionless 
transformation \eqn{\{z^{x}_i := y_i\mathbf{1}_{\{x_i\le x\}}, \,i=1,\cdots,n\}}{z[i]^x={y[i]1(x[i]<=x), i=1,...,n)}}  of 
the observed sample \eqn{\{(x_i,y_i), \,i=1,\cdots,n\}}{(x[i],y[i]), i=1,...,n}, the authors
estimate the conditional endpoint \eqn{\varphi(x)}{varphi(x)} by 
\deqn{\tilde\varphi_{momt}(x) =   z^{x}_{(n-k)}  +  z^{x}_{(n-k)}  M^{(1)}_n \left\{1 + \rho_x  \right\}}{tilde(varphi[momt])(x)= z[n-k]^x + z[n-k]^x M[n]^(1)(1+rho[x])}  
where \eqn{M^{(1)}_n = (1/k)\sum_{i=0}^{k-1}\left(\log  z^x_{(n-i)}- \log   z^x_{(n-k)}\right)}{M[n]^(1)=(1/k)Sum[i=0]^(k-1)(log(z[(n-i)]^x)-log(z[(n-k)]^x))}, 
\eqn{z^{x}_{(1)}\leq \cdots\leq  z^{x}_{(n)}}{z[(1)]^x<=...<=z[(n)]^x} are the ascending order statistics  
corresponding to the transformed sample \eqn{\{z^{x}_i, \,i=1,\cdots,n\}}{z[i]^x, i=1,...,n)}
and \eqn{\rho_x>0}{rho[x]>0}  is referred to as the extreme-value index and has the following interpretation:
when  \eqn{\rho_x>2}{rho[x]>2}, the joint density of data decays smoothly to zero at a speed of power \eqn{\rho_x -2}{rho[x]-2} of the distance from the frontier; 
when \eqn{\rho_x=2}{rho[x]=2},  the density has sudden jumps at the frontier; when \eqn{\rho_x<2}{rho[x]<2}, the density increases toward infinity at a speed of power \eqn{\rho_x -2}{rho[x]-2} 
of the distance from the frontier. 
Most of the contributions to econometric literature on frontier analysis assume that the joint density is strictly positive at its support boundary, or equivalently, \eqn{\rho_x=2}{rho[x]=2} for all \eqn{x}.
When \eqn{\rho_x}{rho[x]} is unknown, Daouia et al. (2010) suggest to use the following two-step estimator: 
First, estimate \eqn{\rho_x}{rho[x]} by the moment estimator \eqn{\tilde\rho_x}{tilde(rho)[x]} implemented in the function \code{\link{rho_momt_pick}} by utilizing the option \code{method="moment"},  
or by the Pickands estimator \eqn{\hat\rho_x}{hat(rho)[x]} by using the option \code{method="pickands"}.
Second, use the estimator \eqn{\tilde\varphi_{momt}(x)}{tilde(varphi[momt])(x)}, as if \eqn{\rho_x}{rho[x]} were known, by substituting the estimated value \eqn{\tilde\rho_x}{tilde(rho)[x]} 
or \eqn{\hat\rho_x}{hat(rho)[x]} in place of \eqn{\rho_x}{rho[x]}.
The \eqn{95\%} confidence interval of \eqn{\varphi(x)}{varphi(x)} derived from the asymptotic normality of \eqn{\tilde\varphi_{momt}(x)}{tilde(varphi[momt])(x)} is given by
\deqn{[\tilde\varphi_{momt}(x)  \pm   1.96 \sqrt{V(\rho_x) / k}  z^{x}_{(n-k)}  M^{(1)}_n (1 + 1/\rho_x) ]}{[tilde(varphi)[momt](x) +- 1.96 sqrt(V(rho[x])/k) z^x[(n-k)]  M^(1)[n](1 + 1/rho[x]) ]}
where \eqn{V(\rho_x) = \rho^2_x  (1+2/\rho_x)^{-1}}{V(rho[x]) = rho^2[x] (1+2/rho[x])^(-1)}.
The sample fraction \eqn{k=k_n(x)}{k=k[n](x)} plays here the role of the smoothing parameter and varies between 1 and \eqn{N_x-1}{N[x]-1}, with \eqn{N_x=\sum_{i=1}^n\mathbf{1}_{\{x_i\le x\}}}{N[x]=Sum[i=1]^n1(x[i]<=x)} 
being the number of observations \eqn{(x_i,y_i)}{(x[i],y[i])} with \eqn{x_i \leq x}{x[i]<=x}. See \code{\link{kopt_momt_pick}} for an automatic data-driven rule for selecting \eqn{k}. 
}

\value{
Returns a numeric vector with the same length as \code{x}.
}


\note{As it is common in extreme-value theory,  good results require a large sample size \eqn{N_x}{N[x]} at each evaluation point \eqn{x}.  
See also the note in \code{\link{kopt_momt_pick}}.}

\references{
Daouia, A., Florens, J.P. and Simar, L. (2010). Frontier Estimation and Extreme Value Theory, \emph{Bernoulli}, 16, 1039-1063.

Dekkers, A.L.M., Einmahl, J.H.J. and L. de Haan (1989), A moment estimator for the index of an extreme-value distribution, \emph{nnals of Statistics}, 17, 1833-1855.
}

\author{
Abdelaati Daouia and Thibault Laurent (converted from Leopold Simar's Matlab code). 
}

\seealso{
\code{\link{dfs_pick}}, \code{\link{kopt_momt_pick}}.
}


\examples{
data("post")
x.post<- seq(post$xinput[100],max(post$xinput), 
 length.out=100) 
# 1. When rho[x] is known and equal to 2, we set:
rho<-2
# To determine the sample fraction k=k[n](x) 
# in tilde(varphi[momt])(x).
best_kn.1<-kopt_momt_pick(post$xinput, post$yprod, 
 x.post, rho=rho)
# To compute the frontier estimates and confidence intervals:  
res.momt.1<-dfs_momt(post$xinput, post$yprod, x.post, 
 rho=rho, k=best_kn.1)
# Representation
plot(yprod~xinput, data=post, xlab="Quantity of labor", 
 ylab="Volume of delivered mail")
lines(x.post, res.momt.1[,1], lty=1, col="cyan")  
lines(x.post, res.momt.1[,2], lty=3, col="magenta")  
lines(x.post, res.momt.1[,3], lty=3, col="magenta")  

\dontrun{
# 2. rho[x] is unknown and estimated by 
# the Pickands estimator tilde(rho[x])
rho_momt<-rho_momt_pick(post$xinput, post$yprod, 
 x.post)
best_kn.2<-kopt_momt_pick(post$xinput, post$yprod,
  x.post, rho=rho_momt)
res.momt.2<-dfs_momt(post$xinput, post$yprod, x.post, 
 rho=rho_momt, k=best_kn.2)  
# 3. rho[x] is unknown independent of x and estimated
# by the (trimmed) mean of tilde(rho[x])
rho_trimmean<-mean(rho_momt, trim=0.00)
best_kn.3<-kopt_momt_pick(post$xinput, post$yprod,
  x.post, rho=rho_trimmean)   
res.momt.3<-dfs_momt(post$xinput, post$yprod, x.post, 
 rho=rho_trimmean, k=best_kn.3)  

# Representation 
plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", 
 ylab="Volume of delivered mail")
lines(x.post, res.momt.2[,1], lty=1, lwd=2, col="cyan")  
lines(x.post, res.momt.2[,2], lty=3, lwd=4, col="magenta")  
lines(x.post, res.momt.2[,3], lty=3, lwd=4, col="magenta")  
plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", 
 ylab="Volume of delivered mail")
lines(x.post, res.momt.3[,1], lty=1, lwd=2, col="cyan")  
lines(x.post, res.momt.3[,2], lty=3, lwd=4, col="magenta")  
lines(x.post, res.momt.3[,3], lty=3, lwd=4, col="magenta") 
}
}

\keyword{nonparametric}