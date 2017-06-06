\name{loc_est_bw}
\alias{loc_est_bw}
\title{
Bandwidth selection for the local linear frontier estimator
}
\description{
Computes the optimal bootstrap bandwidth proposed by Hall and Park (2004) for the local linear frontier estimator. 
}
\usage{
loc_est_bw(xtab, ytab, x, hini, B = 5, method = "u", 
 fix.seed = FALSE, control = list("tm_limit" = 700))
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{hini}{the initial bandwidth at which the local linear estimate will be computed.}
  \item{B}{number of bootstrap replications.}
  \item{method}{a character equal to "u" (unconstrained estimator) or "m" (improved version of the unconstrained estimator).}   
  \item{fix.seed}{a boolean equal to TRUE for fixing the seed (bootstrap sampling).} 
  \item{control}{a list of parameters to the GLPK solver. See *Details* of help(Rglpk_solve_LP).}  
}
\details{
For a detailed description of the bootstrap procedure, see Hall and Park (2004).
}
\note{
The computational burden here is very demanding, so be forewarned. 
}
\value{
Returns the optimal bootstrap bandwidth.
}
\references{
Hall, P. and Park, B.U. (2004). Bandwidth choice for local polynomial estimation of smooth boundaries. \emph{Journal of Multivariate Analysis}, 91, 240-261.
}
\author{
Hohsuk Noh.
}
\seealso{
\code{\link{loc_est}}.
}
\examples{
\dontrun{
data("nuclear")
x.nucl <- seq(min(nuclear$xtab), max(nuclear$xtab), 
 length.out = 101) 
# 1. Unconstrained case 
# Optimal bandwidths over 100 bootstrap replications
system.time(
h.nucl.u <- loc_est_bw(nuclear$xtab, nuclear$ytab, 
 x.nucl, hini = 40, B = 1, method = "u")
)
# result is 79.11877

# 2. Monotonicity constraint
# Optimal bandwidths over 100 bootstrap replications
h.nucl.m <- loc_est_bw(nuclear$xtab, nuclear$ytab, 
 x.nucl, hini = 40, B = 100, method = "m") 
# result is 79.12
}
}

\keyword{optimize}