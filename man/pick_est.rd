\name{pick_est}
\alias{pick_est}
\title{
Local Pickands' frontier estimator
}
\description{
Computes the Pickands type of  estimator introduced by Gijbels and Peng (2000).  
}
\usage{
pick_est(xtab, ytab, x, h, k, type="one-stage")
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{h}{determines the bandwidth at which the estimate will be computed.}
  \item{k}{a numeric vector of the same length as \code{x}, which determines the thresholds at which the Pickands' estimator will be computed.}
  \item{type}{a character equal to "one-stage" or "two-stage".}  
}
\details{
The local Pickands' frontier estimator (option \code{type="one-stage"}), obtained by applying the well-known approach of Dekkers and de Haan (1989) 
in conjunction with the transformed sample of \eqn{z^{xh}_i}{z[i]^(xh)}'s described in the function \code{\link{loc_max}}, is defined as
\deqn{
z^{xh}_{(n-k)} + \left(z^{xh}_{(n-k)}-z^{xh}_{(n-2k)}\right)\{2^{-\log\frac{z^{xh}_{(n-k)}-z^{xh}_{(n-2k)}}{z^{xh}_{(n-2k)}-z^{xh}_{(n-4k)}}/\log 2}-1\}^{-1}.
}{z[n-k]^(xh)+(z[n-k]^(xh)-z[n-2k]^(xh)).(2^(-log((z[n-k]^(xh)-z[n-2k]^(xh))/(z[n-2k]^(xh)-z[n-4k]^(xh)))/log(2))-1)^(-1)}
It is based on three upper order statistics \eqn{z^{xh}_{(n-k)}}{z[n-k]^(xh)}, \eqn{z^{xh}_{(n-2k)}}{z[n-2k]^(xh)}, \eqn{z^{xh}_{(n-4k)}}{z[n-4k]^(xh)}, and depends on \eqn{h} (see \code{\link{loc_max}})  
as well as an intermediate sequence \eqn{k=k(x,n)\to\infty}{k=k(n) tends to infinity} with \eqn{k/n\to 0}{k(n) tends to 0} as \eqn{n\to\infty}{n tends to infinity}. 
The two smoothing parameters \eqn{h} and \eqn{k} have to be fixed in the 4th and 5th arguments of the function.

Also, the user can replace each observation \eqn{y_i}{y[i]} in  the strip of width \eqn{2h} around \eqn{x} by the resulting local Pickands', leaving all observations outside the strip unchanged.
Then, one may apply the DEA estimator (see the function \code{\link{dea_est}}) to the obtained transformed data,
giving the local DEA  estimator  (option  \code{type="two-stage"}). 
}

\value{
Returns a numeric vector with the same length as \code{x}.
}


\references{
Dekkers, A.L.M. and L. de Haan (1989). On the estimation of extreme-value index and large quantiles estimation, \emph{Annals of Statistics}, 17, 1795-1832.

Gijbels, I. and Peng, L. (2000). Estimation of a support curve  via order statistics, \emph{Extremes}, 3,  251-277. 
}

\author{
Abdelaati Daouia and Thibault Laurent.
}

\seealso{
\code{\link{dea_est}}
}


\examples{
\dontrun{
data("green")
plot(log(OUTPUT)~log(COST), data=green)
x <- seq(min(log(green$COST)), max(log(green$COST)), length.out=101)
h=0.5
nx<-unlist(lapply(x,function(y) length(which(abs(log(green$COST)-y)<=h))))
k<-trunc(nx^0.1)
lines(x, pick_est(log(green$COST), log(green$OUTPUT), x, h=h, k=k), lty=1, col="red")
}
}

\keyword{nonparametric}