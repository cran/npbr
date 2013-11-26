\name{loc_est_bw}
\alias{loc_est_bw}
\title{
Bootstrap bandwidth selection for the local linear frontier estimator
}
\description{
Computes the optimal bootstrap bandwidth proposed by Hall and Park (2004) for the local linear frontier estimator 
}
\usage{
loc_est_bw(xtab, ytab, x, hini, B = 5)
}

\arguments{
  \item{xtab}{numeric vectors containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}}
  \item{ytab}{numeric vectors of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed}
  \item{hini}{the initial bandwidth at which the local linear estimate will be computed}
  \item{B}{number of bootstrap replications}
}
\details{
For a detailed description of the bootstrap procedure, see Hall and Park (2004)
}
\value{
Returns the optimal bootstrap bandwidth
}
\references{
Hall, P. and Park, B.U. (2004). Bandwidth choice for local polynomial estimation of smooth boundaries. \emph{Journal of Multivariate Analysis}, 91 (2), 240-261.
}
\author{
Hohsuk Noh
}
\seealso{
\code{\link{loc_est}}
}
\examples{
data("nuclear")
x <- seq(min(nuclear$xtab), max(nuclear$xtab), length.out=1001)
# loc_est_bw(nuclear$xtab, nuclear$ytab, x, hini=40, B=100)
# long computational time 
# returns the value 79.12  
}

\keyword{optimize}