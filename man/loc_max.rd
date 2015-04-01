\name{loc_max}
\alias{loc_max}
\title{
Local maximum frontier estimators
}
\description{
Computes the local constant and local DEA  boundary estimates  proposed by  Gijbels and Peng (2000).   
}
\usage{
loc_max(xtab, ytab, x, h, type="one-stage")
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{h}{determines the bandwidth at which the estimate will be computed.}
  \item{type}{a character equal to "one-stage" or "two-stage".}
}
\details{
When estimating \eqn{\varphi(x)}{varphi(x)}, for a given point \eqn{x\in\R}{x in R}, 
the methodology of Gijbels and Peng consists of considering a strip around \eqn{x} of width \eqn{2h},
where \eqn{h=h_n\to 0}{h=h[n] tend to 0} with \eqn{nh_n\to\infty}{nh[n] tend to infinity}
as \eqn{n\to\infty}{n tend to infinity}, and focusing then on the \eqn{y_i}{y[i]} observations falling into this strip.
More precisely, they consider the transformend variables \eqn{z^{xh}_i = y_i\mathbf{1}_{(|x_i-x|\leq h)}}{z[i]^(xh)=y[i]1[(abs(x[i]-x)<=h)]}, 
\eqn{i=1,\ldots,n}{i=1,...,n}, and the corresponding order statistics \eqn{z^{xh}_{(1)}\le\cdots\le z^{xh}_{(n)}}{z[(1)]^(xh)<=...<=z[(n)]^(xh)}.

The simple maximum \eqn{z^{xh}_{(n)}=\max_{i=1,\ldots,n}z^{xh}_i}{z[(n)]^(xh)=max[i=1,...,n]z[i]^(xh)} defines then the local constant estimator of the
frontier point \eqn{\varphi(x)}{varphi(x)} [option \code{type="one-stage"}].
This opens a way to a  two-stage estimation procedure as follows.
In a first stage, Gijbels and Peng calculate the maximum \eqn{z^{xh}_{(n)}}{z[(n)]^(xh)}.
Then,  they  suggest to replace each observation \eqn{y_i}{y[i]}  in  the strip of width \eqn{2h} around \eqn{x} by  this maximum,  leaving all observations outside the strip unchanged.
More precisely, they define
\eqn{\tilde{y}_i= y_i}{y tilde[i]=y[i]} if \eqn{|x_i-x| > h}{abs(x[i]-x)>h} and \eqn{\tilde{y}_i= z^{xh}_{(n)}}{y tilde[i]=y[i]} if \eqn{|x_i-x| \leq h}{abs(x[i]-x)<=h} either.
Then, they apply the DEA estimator (see the function \code{\link{dea_est}}) to these transformed data \eqn{(x_i,\tilde{y}_i)}{(x[i],tilde y[i])},
giving  the local DEA  estimator  (option  \code{type="two-stage"}). 
An \emph{ad hoc} way of selecting \eqn{h} is by using for instance the function \code{npcdistbw} from the \pkg{np} package (see Daouia et al. (2015) for details).
}

\value{
Returns a numeric vector with the same length as \code{x}.
}


\references{
Daouia, A., Laurent, T. and Noh, H. (2015). \pkg{npbr}: A Package for Nonparametric Boundary Regression in R. 

Gijbels, I. and Peng, L. (2000). Estimation of a support curve  via order statistics, \emph{Extremes}, 3,  251--277. 
}

\author{
Abdelaati Daouia and Thibault Laurent.
}

\seealso{
\code{\link{dea_est}}
}


\examples{
data("green")
x.green <- seq(min(log(green$COST)), max(log(green$COST)), 
 length.out=101)
# Local maximum frontier estimates
# a. Local constant estimator
loc_max_1stage<-loc_max(log(green$COST), log(green$OUTPUT), 
 x.green, h=0.5, type="one-stage")
# b. Local DEA estimator
loc_max_2stage<-loc_max(log(green$COST), log(green$OUTPUT), 
 x.green, h=0.5, type="two-stage")  
# Representation 
plot(log(OUTPUT)~log(COST), data=green)
lines(x.green, loc_max_1stage, lty=1, col="magenta")
lines(x.green, loc_max_2stage, lty=2, col="cyan")
legend("topleft",legend=c("one-stage", "two-stage"), 
 col=c("magenta","cyan"), lty=c(1,2))
}

\keyword{nonparametric}