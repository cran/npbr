\name{kopt_momt_pick}
\alias{kopt_momt_pick}
\title{
Optimal \eqn{k} in moment and Pickands frontier estimators
}
\description{
This function gives the optimal sample fraction k in the moment and Pickands type of estimators introduced by Daouia, Florens and Simar (2010).    
}
\usage{
kopt_momt_pick(xtab, ytab, x, rho, method="moment", wind.coef=0.1)
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{rho}{a numeric vector of the same length as \code{x} or a scalar, which determines the values of rho.}
  \item{method}{a character equal to "moment" or "pickands".}  
  \item{wind.coef}{a scalar coefficient to be selected in the interval (0,1].}    
}

\details{
This function is an implementation of an experimental method by Daouia et al. (2010) for the 
automated threshold selection (choice of \eqn{k=k_n(x)}{k=k[n](x)}) for the moment frontier estimator 
\eqn{\tilde\varphi_{momt}(x)}{tilde(varphi)[momt](x)}  [see \code{\link{dfs_momt}}] in case  \code{method="moment"} and for the Pickands 
frontier estimator \eqn{\hat\varphi_{pick}(x)}{hat(varphi)[pick](x)} [see \code{\link{dfs_pick}}] in case \code{method="pickands"}.
The idea is to select first (for each \eqn{x}) a grid of values for the sample fraction \eqn{k_n(x)}{k[n](x)} given by \eqn{k = 1, \cdots, [\sqrt{N_x}]}{k=1,...,[sqrt(N[x])]}, 
where \eqn{[\sqrt{N_x}]}{[sqrt(N[x])]} stands for the integer part of \eqn{\sqrt{N_x}}{sqrt(N[x])} with \eqn{N_x=\sum_{i=1}^n1_{\{x_i\le x\}}}{N[x]=sum[i=1]^n 1(x[i]<=x)}, 
and then select the \eqn{k} where the variation of the results is the smallest. To achieve this here, 
Daouia et al. (2010) compute the standard deviations of \eqn{\tilde\varphi_{momt}(x)}{tilde(varphi)[momt](x)} [option \code{method="moment"}] or \eqn{\hat\varphi_{pick}(x)}{hat(varphi)[pick](x)} [option \code{method="pickands"}]
over a ``window'' of size \eqn{\max(3, [ wind.coef  \times \sqrt{N_x} /2])}{max(3,[wind.coeff x sqrt(N[x])/2]}, where the coefficient \code{wind.coef} should be selected in the interval \eqn{(0,1]} in such 
a way to avoid numerical instabilities. 
The default option \code{wind.coef=0.1} corresponds to having a window large enough to cover around \eqn{10\%} of the possible values of \eqn{k} in the selected range of values for \eqn{k_n(x)}{k[n](x)}. 
The value of \eqn{k} where the standard deviation is minimal defines the desired sample fraction \eqn{k_n(x)}{k[n](x)}. 
}

\note{
In order to choose a raisonable estimate  \eqn{\tilde\varphi_{momt}(x)=\tilde\varphi_{momt}(x,k)}{tilde(varphi)[momt](x)=tilde(varphi)[momt](x,k)} [see \code{\link{dfs_momt}}] and 
\eqn{\hat\varphi_{pick}(x)=\hat\varphi_{pick}(x,k)}{hat(varphi)[pick](x)=hat(varphi)[pick](x,k)} [see \code{\link{dfs_pick}}] of the frontier function  \eqn{\varphi(x)}{varphi(x)}, 
for each fixed \eqn{x}, one can construct the plot of the estimator of interest, consisting of the points \eqn{\{(k,\tilde\varphi_{momt}(x,k))\}_k}{{(k,tilde(varphi)[momt](x,k))}[k]} or
\eqn{\{(k,\hat\varphi_{pick}(x,k))\}_k}{{(k,hat(varphi)[pick](x,k))}[k]}, and select a value of the estimate at which the obtained graph looks stable. This is this kind of idea
which guides the proposed automatic data-driven rule for a chosen grid of values of \eqn{x}. The main difficulty with such a method is that the plots of
\eqn{\tilde\varphi_{momt}(x,k)}{tilde(varphi)[momt](x,k)} or \eqn{\hat\varphi_{pick}(x,k)}{hat(varphi)[pick](x,k)} as functions of \eqn{k}, for each \eqn{x}, may be so unstable that reasonable values of
\eqn{k} [which would correspond to the true value of \eqn{\varphi(x)}{varphi(x)}] may be hidden in the graphs. In result, the obtained frontier estimator may exhibit considerable volatility as
a function of \eqn{x}. One way to avoid such instabilities is by tuning the choice of the parameter \code{wind.coef} in the interval (0,1]. Note that the default value is \code{wind.coef=0.1}.
The user can also improve appreciably the estimation of  \eqn{\varphi(x)}{varphi(x)} by refining the estimation of the extreme-value index \eqn{\rho_x}{rho[x]} (see \code{\link{rho_momt_pick}} for details).  
}

\value{
Returns a numeric vector with the same length as \code{x}.
}


\references{
Daouia, A., Florens, J.P. and Simar, L. (2010). Frontier Estimation and Extreme Value Theory, \emph{Bernoulli}, 16, 1039-1063.

Dekkers, A.L.M., Einmahl, J.H.J. and L. de Haan (1989), A moment estimator for the index of an extreme-value distribution, \emph{Annals of Statistics}, 17, 1833-1855.
}

\author{
Abdelaati Daouia and Thibault Laurent (converted from Leopold Simar's Matlab code). 
}

\seealso{
\code{\link{dfs_momt}}, \code{\link{dfs_pick}}.
}


\examples{
data("post")
x.post<- seq(post$xinput[100],max(post$xinput), 
 length.out=100) 
# When rho[x] is known and equal to 2, we set:
rho<-2
# a. Optimal k in Pickands frontier estimators
best_kn.pick<-kopt_momt_pick(post$xinput, post$yprod, 
 x.post, method="pickands", rho=rho)
# b. Optimal k in moment frontier estimators
\dontrun{
best_kn.momt<-kopt_momt_pick(post$xinput, post$yprod, 
 x.post, rho=rho)
}
}

\keyword{nonparametric}