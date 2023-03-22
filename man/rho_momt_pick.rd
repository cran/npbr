\name{rho_momt_pick}
\alias{rho_momt_pick}
\title{
Optimal rho for moment and Pickands frontier estimator
}
\description{
This function gives the optimal rho involved in the moment and Pickands estimators of Daouia, Florens and Simar (2010).    
}
\usage{
rho_momt_pick(xtab, ytab, x, method="moment", lrho=1, urho=Inf)
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{method}{a character equal to "moment" or "pickands".}
  \item{lrho}{a scalar, minimum rho threshold value.}  
  \item{urho}{a scalar, maximum rho threshold value.}
}
\details{
This function computes the moment and Pickands estimates of the extreme-value index 
\eqn{\rho_x}{rho[x]} involved in the frontier estimators \eqn{\tilde\varphi_{momt}(x)}{tilde(varphi)[momt](x)}  [see \code{\link{dfs_momt}}] and 
\eqn{\hat\varphi_{pick}(x)}{hat(varphi)[pick](x)} [see \code{\link{dfs_pick}}].
In case \code{method="moment"}, the estimator of \eqn{\rho_x}{rho[x]} defined as
\deqn{\tilde{\rho}_x = -\left(M^{(1)}_n + 1 -\frac{1}{2}\left[1-(M^{(1)}_n)^2/M^{(2)}_n\right]^{-1}\right)^{-1}}{tilde(rho)[x] = -(M^(1)[n] + 1 - [1-(M^(1)[n])^2/M^(2)[n]]^(-1)/2)^(-1)}
is based on the moments \eqn{M^{(j)}_n = (1/k)\sum_{i=0}^{k-1}\left(\log  z^x_{(n-i)}- \log   z^x_{(n-k)}\right)^j}{M[n]^(j)=(1/k)Sum[i=0]^(k-1)(log(z[n-i]^x)-log(z[n-k]^x))^j}   
for \eqn{j=1,2}, with \eqn{z^{x}_{(1)}\leq \cdots\leq  z^{x}_{(n)}}{z[(1)]^x<=...<=z[(n)]^x} are the ascending order statistics  
corresponding to the transformed sample \eqn{\{z^{x}_i := y_i\mathbf{1}_{\{x_i\le x\}}, \,i=1,\cdots,n\}}{z[i]^x={y[i]1(x[i]<=x), i=1,...,n)}}
In case \code{method="pickands"}, the estimator of \eqn{\rho_x}{rho[x]} is given by
\deqn{\hat{\rho}_x = - \log 2/\log\{(z^x_{(n-k+1)} - z^x_{(n-2k+1)})/(z^x_{(n-2k+1)} - z^x_{(n-4k+1)})\}.}{hat(rho)[x] = - log(2)/log((z^x[(n-k+1)] - z^x[(n-2k+1)])/(z^x[(n-2k+1)] - z^x[(n-4k+1)])).}
To select the threshold \eqn{k=k_n(x)}{k=k[n](x)} in \eqn{\tilde{\rho}_x}{tilde(rho)[x]} and \eqn{\hat{\rho}_x}{hat(rho)[x]}, Daouia et al. (2010) have suggested to use the following data driven method for each 
\eqn{x}: They first select a grid of values for \eqn{k=k_n(x)}{k=k[n](x)}.
For the Pickands estimator \eqn{\hat{\rho}_x}{hat(rho)[x]}, they choose \eqn{k_n(x) = [N_x /4] - k + 1}{k[n](x) = [N[x] /4] - k + 1}, where \eqn{k} is an integer varying between 1 
and the integer part \eqn{[N_x/4]}{[N[x]/4]} of \eqn{N_x/4}{N[x]/4}, with \eqn{N_x=\sum_{i=1}^n1_{\{x_i\le x\}}}{N[x]=sum[i=1]^n 1[x[i] <= x]}.
For the moment estimator \eqn{\tilde{\rho}_x}{tilde(rho)[x]}, they choose \eqn{k_n(x) = N_x - k}{k[n](x) = N[x] - k}, where \eqn{k} is an integer varying between 1 and \eqn{N_x -1}{N[x]-1}.
Then, they evaluate the estimator \eqn{\hat{\rho}_x(k)}{hat(rho)[x](k)}  (respectively, \eqn{\tilde{\rho}_x(k)}{tilde(rho)[x](k)}) and select the k where the variation of the results is the smallest. 
They achieve this by computing the standard deviation of \eqn{\hat{\rho}_x(k)}{hat(rho)[x](k)} (respectively, \eqn{\tilde{\rho}_x(k)}{tilde(rho)[x](k)}) over a ``window'' of 
\eqn{\max([\sqrt{N_x /4}],3)}{max(sqrt(N[x]/4),3)} (respectively, \eqn{\max([\sqrt{N_x-1}],3)}{max(sqrt(N[x]-1),3)}) 
successive values of \eqn{k}. The value of \eqn{k} where this standard deviation is minimal defines the value of \eqn{k_n(x)}{k[n](x)}.
The user can also appreciably improve the estimation of \eqn{\rho_x}{rho[x]} and \eqn{\varphi(x)}{varphi(x)} itself by tuning the choice of the lower limit (default option \code{lrho=1}) 
and upper limit (default option \code{urho=Inf}).
}

\note{
In order to choose a raisonable estimate \eqn{\tilde\rho_x=\tilde\rho_x(k)}{tilde(rho)[x]=tilde(rho)[x](k)} and 
\eqn{\hat\rho_x=\hat\rho_x(k)}{hat(rho)[x]=hat(rho)[x](k)} of the extreme-value index \eqn{\rho_x}{rho[x]}, 
for each fixed \eqn{x}, one can construct the plot of the estimator of interest, consisting of the points \eqn{\{(k,\tilde\rho_x(k))\}_k}{{(k,tilde(rho)[x](k))}[k]} or
\eqn{\{(k,\hat\rho_x(k))\}_k}{{(k,hat(rho)[x](k))}[k]}, and select a value of the estimate at which the obtained graph looks stable. This is this kind of idea
which guides the propoed automatic data-driven rule for a chosen grid of values of \eqn{x}. The main difficulty with such a method is that the plots of
\eqn{\tilde\rho_x(k)}{tilde(rho)[x](k)} or \eqn{\hat\rho_x(k)}{hat(rho)[x](k)} as functions of \eqn{k}, for each \eqn{x}, may be so unstable that reasonable values of
\eqn{k} [which would correspond to the true value of \eqn{\rho_x}{rho[x]}] may be hidden in the graphs. In results, the obtained extreme-value index estimator and the frontier estimator itself may 
exhibits considerable volatility as functions of \eqn{x}. The user can appreciably improve the estimation of \eqn{\rho_x}{rho[x]} and \eqn{\varphi(x)}{varphi(x)} 
by tuning the choice of the lower limit (default option \code{lrho=1}) and upper limit (default option \code{urho=Inf}). 
}

\value{
Returns a numeric vector with the same length as \code{x}.
}

\references{
Daouia, A., Florens, J.P. and Simar, L. (2010). Frontier Estimation and Extreme Value Theory, \emph{Bernoulli}, 16, 1039-1063.

Dekkers, A.L.M., Einmahl, J.H.J. and L. de Haan (1989), A moment estimator for the index of an extreme-value distribution, \emph{The Annals of Statistics}, 17(4), 1833-1855.
}

\author{
Abdelaati Daouia and Thibault Laurent (codes converted from Matlab's Leopold Simar code). 
}

\seealso{
\code{\link{dfs_momt}}, \code{\link{dfs_pick}}
}

\examples{
data("post")
x.post<- seq(post$xinput[100],max(post$xinput), 
 length.out=100) 
\dontrun{
# a. Optimal rho for Pickands frontier estimator
rho_pick<-rho_momt_pick(post$xinput, post$yprod, 
 x.post, method="pickands")
# b. Optimal rho for moment frontier estimator
rho_momt<-rho_momt_pick(post$xinput, post$yprod, 
 x.post, method="moment")
}
}

\keyword{nonparametric}