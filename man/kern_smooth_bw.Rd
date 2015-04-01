\name{kern_smooth_bw}
\alias{kern_smooth_bw}
\title{
Bandwidth selection for kernel smoothing frontier estimators
}
\description{
The function \code{kern_smooth_bw} provides two bandwidth selection methods. One is the least squares cross-validation developed by Parmeter and Racine (2013). The other is the BIC developed in Noh (2014).    
}
\usage{
kern_smooth_bw(xtab, ytab, method="u", technique="noh", bw_method="bic")
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{method}{a character equal to "u" (unconstrained estimator), "m" (under the monotonicity constraint) or "mc" (under simultaneous monotonicity and concavity constraints).}  
  \item{technique}{which estimation technique to use: "Noh" specifies the use of the method in Noh (2014), while "pr" is for the method in Parameter and Racine (2013).}  
  \item{bw_method}{which bandwidth selection method to use: "cv" returns the bandwidth that minimizes the least squares cross-validation criterion, and "bic" returns the bandwidth minimizing the BIC.}
}
\details{
As with any smoothed techniques, the bandwidth selection is critical to the quality of the frontier estimator. Parmeter and Racine (2013)'s recommendation is to use the least squares cross-validation method implemented with \code{bw\_method="cv"} in the function \code{kern\_smooth\_bw}.
Instead, Noh (2014) proposed to select the bandwidth which minimizes the following criterion: 
\deqn{ 
BIC(h) = \log \left( \sum_{i=1}^n (\hat \varphi(x_i|\hat p(h))-y_i)\right)+\frac {\log n \cdot tr(S(h))}{2n},}
where \eqn{\hat p(h)} is the chosen weight vector associated to the bandwidth \eqn{h}, and \eqn{tr(S(h))} is the trace of the smoothing matrix 
\deqn{
S(h) = \left( \begin{array}{ccc} A_1(x_1) & \cdots & A_n(x_1)
\\ \vdots & \ddots& \vdots
\\ A_1(x_n) & \cdots & A_n(x_n) \end{array} \right).
}
The function \code{kern\_smooth\_bw} computes the optimal bandwidth from this criterion with option \code{bw\_method="bic"}.
}

\value{
Returns an optimal bandwidth depending on the specified selection method. 
}

\references{
Noh, H. (2014). Frontier estimation using kernel smoothing estimators with data transformation. \emph{Journal of the Korean Statistical Society}, 43, 503-512.

Parmeter, C.F. and Racine, J.S. (2013). Smooth constrained frontier analysis in \emph{Recent Advances and Future Directions in Causality, Prediction, and Specification Analysis}, Springer-Verlag, New York, 463-488.
}

\author{Hohsuk Noh
}

\seealso{
\code{\link{kern_smooth}}.
}

\examples{
data("green")
x.green <- seq(min(log(green$COST)), max(log(green$COST)),length.out=101)
options(np.tree=TRUE,crs.messages=FALSE,np.messages=FALSE)
h.pr.green.m<-kern_smooth_bw(log(green$COST),log(green$OUTPUT), method="m", 
 technique="pr", bw_method="cv")
h.noh.green.m<-kern_smooth_bw(log(green$COST),log(green$OUTPUT), method="m", 
 technique="noh", bw_method="bic")
y.pr.green.m<-kern_smooth(log(green$COST),log(green$OUTPUT), x.green, 
 h=h.pr.green.m, method="m", technique="pr")
y.noh.green.m<-kern_smooth(log(green$COST),log(green$OUTPUT), x.green, 
 h=h.noh.green.m, method="m", technique="noh")
plot(log(OUTPUT)~log(COST), data=green, xlab="log(COST)",ylab="log(OUTPUT)") 
lines(x.green, y.pr.green.m, lwd=4, lty=3, col="red") 
lines(x.green, y.noh.green.m, lwd=4, lty=3, col="blue")  
legend("topleft", col=c("blue","red"),lty=3, legend=c("noh","pr"), 
 lwd=4, cex=0.8)
}
\keyword{optimize}