\name{dfs_pick}
\alias{dfs_pick}
\title{
Pickands frontier estimator
}
\description{
This function is an implementation of the Pickands type estimator developed by Daouia, Florens and Simar (2010).    
}
\usage{dfs_pick(xtab, ytab, x, k, rho, ci=TRUE)
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}  
  \item{k}{a numeric vector of the same length as \code{x} or a scalar, which determines the thresholds at which the Pickands estimator will be computed.}
  \item{rho}{a numeric vector of the same length as \code{x} or a scalar, which determines the values of rho.}
  \item{ci}{a boolean, TRUE for computing the confidence interval.}    
}
\details{
Built on the ideas of Dekkers and de Haan (1989), Daouia et al. (2010) propose to estimate the frontier point \eqn{\varphi(x)}{varphi(x)} by
\deqn{\hat\varphi_{pick}(x) =   \frac{z^x_{(n-k+1)} - z^x_{(n-2k+1)}}{2^{1/\rho_x} - 1} + z^x_{(n-k+1)}}{hat(varphi)[pick](x) =   (z^x[(n-k+1)] - z^x[(n-2k+1)])/(2^(1/rho[x] - 1) + z^x_[(n-k+1)]}
from the transformed data \eqn{\{z^{x}_i, \,i=1,\cdots,n\}}{z[i]^x, i=1,...,n} described in \code{\link{dfs_momt}}, where \eqn{\rho_x>0}{rho[x]>0} is the same tail-index as in \code{\link{dfs_momt}}.
If \eqn{\rho_x}{rho[x]}  is known (typically equal to 2 if the joint density of data is believed to have sudden jumps at the frontier), then one can use the estimator \eqn{\hat\varphi_{pick}(x)}{hat(varphi)[pick](x)}
in conjunction with the function \code{\link{kopt_momt_pick}} which implements an automatic data-driven method for selecting the threshold \eqn{k}. 
In contrast, if \eqn{\rho_x}{rho[x]}  is unknown, one could consider using the following two-step estimator:
First, estimate \eqn{\rho_x}{rho[x]} by the Pickands estimator \eqn{\hat\rho_x}{hat(rho)[x]}  implemented in the function \code{\link{rho_momt_pick}} by using the option \code{method="pickands"}, 
or by the moment estimator \eqn{\tilde\rho_x}{tilde(rho)[x]} by utilizing the option \code{method="moment"}. 
Second, use the estimator \eqn{\hat\varphi_{pick}(x)}{hat(varphi)[pick](x)}, as if \eqn{\rho_x}{rho[x]} were known, by substituting the estimated value \eqn{\hat\rho_x}{hat(rho)[x]} or \eqn{\tilde\rho_x}{tilde(rho)[x]} 
in place of \eqn{\rho_x}{rho[x]}.
The pointwise \eqn{95\%} confidence interval of the frontier function obtained from the asymptotic normality of \eqn{\hat\varphi_{pick}(x)}{hat(varphi)[pick](x)} is given by
\deqn{[\hat\varphi_{pick}(x)  \pm   1.96 \sqrt{v(\rho_x) / (2 k)}    ( z^x_{(n-k+1)} - z^x_{(n-2k+1)})]}{[hat(varphi)[pick](x) +- 1.96 sqrt(v(rho[x])/2k) (z^x[(n-k+1)]  - z^x[(n-2k+1)])]}
where \eqn{v(\rho_x) =\rho^{-2}_x  2^{-2/\rho_x}/(2^{-1/\rho_x} -1)^4}{v(rho[x]) =rho^(-2)[x]  2^(-2/rho[x])/(2^(-1/rho[x]) -1)^4}.
Finally, to select the threshold \eqn{k=k_n(x)}{k=k[n](x)}, one could use the automatic data-driven method of Daouia et al. (2010) implemented in the function 
\code{\link{kopt_momt_pick}} (option \code{method="pickands"}).
}

\value{
Returns a numeric vector with the same length as \code{x}.
}


\note{As it is common in extreme-value theory,  good results require a large sample size \eqn{N_x}{N[x]} at each evaluation point \eqn{x}. See also the note in \code{\link{kopt_momt_pick}}.}


\references{
Daouia, A., Florens, J.P. and Simar, L. (2010). Frontier Estimation and Extreme Value Theory, \emph{Bernoulli}, 16, 1039-1063.

Dekkers, A.L.M., Einmahl, J.H.J. and L. de Haan (1989), A moment estimator for the index of an extreme-value distribution, \emph{Annals of Statistics}, 17, 1833-1855.
}

\author{
Abdelaati Daouia and Thibault Laurent (converted from Leopold Simar's Matlab code). 
}

\seealso{
\code{\link{dfs_momt}}, \code{\link{kopt_momt_pick}}.
}


\examples{
data("post")
x.post<- seq(post$xinput[100],max(post$xinput), 
 length.out=100) 
# 1. When rho[x] is known and equal to 2, we set:
rho<-2
# To determine the sample fraction k=k[n](x) 
# in hat(varphi[pick])(x).
best_kn.1<-kopt_momt_pick(post$xinput, post$yprod, 
 x.post, method="pickands", rho=rho)
# To compute the frontier estimates and confidence intervals:  
res.pick.1<-dfs_pick(post$xinput, post$yprod, x.post, 
 rho=rho, k=best_kn.1)
# Representation
plot(yprod~xinput, data=post, xlab="Quantity of labor", 
 ylab="Volume of delivered mail")
lines(x.post, res.pick.1[,1], lty=1, col="cyan")  
lines(x.post, res.pick.1[,2], lty=3, col="magenta")  
lines(x.post, res.pick.1[,3], lty=3, col="magenta")  

\dontrun{
# 2. rho[x] is unknown and estimated by 
# the Pickands estimator hat(rho[x])
rho_pick<-rho_momt_pick(post$xinput, post$yprod, 
 x.post, method="pickands")
best_kn.2<-kopt_momt_pick(post$xinput, post$yprod,
  x.post, method="pickands", rho=rho_pick)
res.pick.2<-dfs_pick(post$xinput, post$yprod, x.post, 
 rho=rho_pick, k=best_kn.2)  
# 3. rho[x] is unknown independent of x and estimated
# by the (trimmed) mean of hat(rho[x])
rho_trimmean<-mean(rho_pick, trim=0.00)
best_kn.3<-kopt_momt_pick(post$xinput, post$yprod,
  x.post, rho=rho_trimmean, method="pickands")   
res.pick.3<-dfs_pick(post$xinput, post$yprod, x.post, 
 rho=rho_trimmean, k=best_kn.3)  

# Representation 
plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", 
 ylab="Volume of delivered mail")
lines(x.post, res.pick.2[,1], lty=1, lwd=2, col="cyan")  
lines(x.post, res.pick.2[,2], lty=3, lwd=4, col="magenta")  
lines(x.post, res.pick.2[,3], lty=3, lwd=4, col="magenta")  
plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", 
 ylab="Volume of delivered mail")
lines(x.post, res.pick.3[,1], lty=1, lwd=2, col="cyan")  
lines(x.post, res.pick.3[,2], lty=3, lwd=4, col="magenta")  
lines(x.post, res.pick.3[,3], lty=3, lwd=4, col="magenta") 
}
}

\keyword{nonparametric}