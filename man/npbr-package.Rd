\name{npbr-package}
\alias{npbr-package}
\alias{npbr}
\docType{package}
\title{Nonparametric boundary regression}
\description{
This package provides a variety of nonparemtric regression methods and allows for comparison among these methods via Monte Carlo experiments. The package provides also illustrations with real datasets.}
\details{
Suppose that we have \eqn{n} pairs of observations
\eqn{(x_i,y_i),~i=1,\ldots,n}{(xi,yi), i=1,...,n},
from a bivariate distribution with a density \eqn{f(x,y)} in \eqn{R^2}{R2}. The support \eqn{\Psi}{Psi} of \eqn{f} is assumed to be of the form
  \deqn{
 \Psi = \{ (x,y) | y \leq \varphi(x) \}   \supseteq   \{ (x,y) | f(x,y) > 0 \} }{
  Psi = {(x,y), y <= varphi(x)}  supseteq   {(x,y) s.a. f(x,y) > 0 }
  }
  \deqn{\{ (x,y) | y > \varphi(x) \}   \subseteq \{ (x,y) | f(x,y) = 0 \},
  }{{(x,y), y > varphi(x)} supseteq  {(x,y) s.a. f(x,y) = 0 }}  
  
where \eqn{\varphi}{varphi} is a monotone increasing
and/or concave function whose graph corresponds to the locus of the curve above which the density \eqn{f} is zero.
We consider the estimation of the frontier function \eqn{\varphi}{varphi} based on the sample \eqn{\{ (x_i,y_i),~i=1,\ldots,n\}}{{(x_i,y_i), i=1,...,n}}
in the general setting where the density \eqn{f} may have sudden jumps at the frontier, decay to zero or rise up to infinity as it approaches its support boundary.

The package provides functions for the best known nonparametric estimation procedures.
The selected methods can be divided into a number of different categories: empirical, smooth, unconstrained and restricted estimates.
The package provides some real data sets as well.
}
\author{
Abdelaati Daouia <Abdelaati.Daouia@tse-fr.eu>, Thibault Laurent <thibault.laurent@univ-tlse1.fr>, Hohsuk Noh <word5810@gmail.com>

Maintainer: Thibault Laurent <thibault.laurent@univ-tlse1.fr>
}
\references{
Daouia, A., Noh, H. and Park, B.U. (2013). Data Envelope fitting with constrained polynomial splines. \emph{TSE Working Paper}, \url{http://www.tse-fr.eu/images/doc/wp/etrie/wp_tse_449.pdf}.
}

\keyword{nonparametric}
\keyword{optimize}

\examples{
data("green")
plot(log(OUTPUT)~log(COST), data=green)
x <- seq(min(log(green$COST)), max(log(green$COST)), length.out=1001)
# under the separate monotonicity constraint with the knots mesh being
# obtained via AIC and BIC criteria
lines(x, quad_spline_est(log(green$COST), log(green$OUTPUT), x, kn=9, cv=0),
lty=1, col="red")   
lines(x, quad_spline_est(log(green$COST), log(green$OUTPUT), x, kn=6, cv=0),
lty=1, col="blue")
# under both monotonicity and concavity constraints
# with the same optimal number of knots obtained using AIC and BIC criteria
lines(x, quad_spline_est(log(green$COST), log(green$OUTPUT), x, kn=1, cv=1), 
lty=2)   
legend("topleft",legend=c("monotonicity (kn=9)", "monotonicity (kn=6)",
"monotonicity + concavity (kn=1)"), col=c("red","blue","black"), lty=c(1,1,2))
}
