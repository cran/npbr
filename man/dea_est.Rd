\name{dea_est}
\alias{dea_est}
\title{
DEA, FDH and linearized FDH estimators.  
}
\description{
The function implements the empirical FDH (free disposal hull), LFDH (linearized FDH) and DEA (data envelopment analysis) frontier estimators.  
}
\usage{
dea_est(xtab, ytab, x, type = "dea")
}

\arguments{
  \item{xtab}{numeric vectors containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}}
  \item{ytab}{numeric vectors of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed}
  \item{type}{a character equal to "dea", "fdh" or "lfdh"}
}
\details{
There are mainly two usual frontier estimation methods for preserving monotonicity: the free disposal hull (FDH) introduced by Deprins \emph{et al.} (1984)  
and the data envelopment analysis (DEA) initiated by Farrell (1957). The FDH boundary is the lowest ``stair-case'' monotone curve covering all the data points
\deqn{\varphi_n(x):=\max\{y_i,\,i:x_i\leq x\}.}{varphi[n](x):=max{yi, i:xi <= x}.}
An improved version of this estimator, referred to as the linearized FDH (LFDH), is obtained by drawing the polygonal line smoothing the staircase FDH curve. 
It has been considered in Hall and Park (2002) and Jeong and Simar (2006).
When the joint support of data is in addition convex, the DEA estimator is defined as the least concave majorant of the FDH frontier
(see also Gijbels \emph{et al.} (1999)).

Employs the function \code{DEA} from the package \pkg{Benchmarking}.
}

\value{
Returns a numeric vector with the same length as \code{x}
}
\references{
Bogetoft, P. and Otto, L. (2011), \emph{Benchmarking with DEA, SFA and R}, Springer-Verlag   

Deprins, D., Simar, L. and H. Tulkens (1984).  Measuring labor
  efficiency in post offices, in: M. Marchand, P. Pestieau and H. Tulkens,
  (Eds), \emph{The performance of Public Enterprises: Concepts and Measurements}. North-Holland, Amsterdam, pp. 243--267. 

Farrell, M.J. (1957). The measurement of productive efficiency.
  \emph{Journal of the Royal Statistical Society, Series A}, 120, 253--281. 

Gijbels, I., Mammen, E., Park, B.U. and Simar, L. (1999). On estimation of monotone and concave frontier functions, \emph{Journal of American Statistical Association}, 94, 220--228.   

Hall, P. and Park, B.U. (2002). New methods for bias correction at endpoints and boundaries, \emph{Annals of Statistics}, 30, 1460-1479.  

Jeong, S.-O. and Simar, L. (2006). Linearly interpolated FDH efficiency score for nonconvex frontiers,
\emph{Journal of Multivariate Analysis}, 97, 2141--2161.
}
\author{
Hohsuk Noh
}

\seealso{
\code{\link{quad_spline_est}}
}

\examples{
data("green")
plot(OUTPUT~COST, data=green)
x <- seq(min(green$COST), max(green$COST), length.out=1001)
# dea 
lines(x, dea_est(green$COST, green$OUTPUT, x, type="dea"),
lty=1, col="red")  
# fdh 
lines(x, dea_est(green$COST, green$OUTPUT, x, type="fdh"),
lty=2, col="blue")
# lfdh
lines(x, dea_est(green$COST, green$OUTPUT, x, type="lfdh"), 
lty=3, col="green")   
legend("topleft",legend=c("dea","fdh","lfdh"), col=c("red","blue","green"), 
lty=1:3)
}
\keyword{nonparametric}
\keyword{optimize}
