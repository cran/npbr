\name{kern_smooth}
\alias{kern_smooth}
\title{
Frontier estimation via kernel smoothing  
}
\description{
The function \code{kern_smooth} implements two frontier estimators based on kernel smoothing techniques. One is from Noh (2014) and the other is from Parmeter and Racine (2013).  
}
\usage{
kern_smooth(xtab, ytab, x, h, method="u", technique="noh", 
 control = list("tm_limit" = 700))
}

\arguments{
  \item{xtab}{a numeric vector containing the observed inputs  \eqn{x_1,\ldots,x_n}{x1,...,xn}.}
  \item{ytab}{a numeric vector of the same length as \code{xtab} containing the observed outputs \eqn{y_1,\ldots,y_n}{y1,...,yn}.}
  \item{x}{a numeric vector of evaluation points in which the estimator is to be computed.}
  \item{h}{determines the bandwidth at which the smoothed kernel estimate will be computed.}
  \item{method}{a character equal to "u" (unconstrained estimator), "m" (under the monotonicity constraint) or "mc" (under simultaneous monotonicity and concavity constraints).}  
  \item{technique}{which estimation method to use. "Noh"" specifies the use of the method in Noh (2014) and "pr" is for the method in Parameter and Racine (2013).}
  \item{control}{a list of parameters to the GLPK solver. See *Details* of help(Rglpk_solve_LP).}  
  }  

\details{
To estimate the frontier function, Parameter and Racine (2013) considered the following generalization of linear regression smoothers 
\deqn{\hat \varphi(x|p) = \sum_{i=1}^n p_i A_i(x)y_i,} where \eqn{A_i(x)} is the kernel weight function of \eqn{x} for the \eqn{i}th data depending on \eqn{x_i}'s and the sort of linear smoothers. For example, the Nadaraya-Watson kernel weights are \eqn{A_i(x) = K_i(x)/(\sum_{j=1}^n K_j(x)),} where \eqn{K_i(x) = h^{-1} K\{ (x-x_i)/h\}}, with the kernel function \eqn{K} being a bounded and symmetric probability density, and \eqn{h} is a bandwidth. Then, the weight vector \eqn{p=(p_1,\ldots,p_n)^T} is chosen to minimize the distance \eqn{D(p)=(p-p_u)^T(p-p_u)} subject to the envelopment constraints and the choice of the shape constraints, where \eqn{p_u} is an \eqn{n}-dimensional vector with all elements being one. The envelopement and shape constraints are 
\deqn{
\begin{array}{rcl}
\hat \varphi(x_i|p) - y_i &=& \sum_{i=1}^n p_i A_i(x_i)y_i - y_i \geq 0,~i=1,\ldots,n;~~~{\sf (envelopment~constraints)} 
\\ \hat \varphi^{(1)}(x|p) &=& \sum_{i=1}^n p_i A_i^{(1)}(x)y_i \geq 0,~x \in \mathcal{M};~~~{\sf (monotonocity~constraints)} 
\\  \hat \varphi^{(2)}(x|p) &=& \sum_{i=1}^n p_i A_i^{(2)}(x)y_i \leq 0,~x \in \mathcal{C},~~~{\sf (concavity~constraints)} 
\end{array}}
where \eqn{\hat \varphi^{(s)}(x|p) = \sum_{i=1}^n p_i A_i^{(s)}(x) y_i} is the \eqn{s}th derivative of \eqn{\hat \varphi(x|p)}, with \eqn{\mathcal{M}} and \eqn{\mathcal{C}} being the collections of points where monotonicity and concavity are imposed, respectively. In our implementation of the estimator, we simply take the entire dataset \eqn{ \{(x_i,y_i),~i=1,\ldots,n\}} to be \eqn{\mathcal{M}} and \eqn{\mathcal{C}} and, in case of small samples, we augment the sample points by an equispaced grid of length 201 over the observed support \eqn{[\min_i x_i,\max_i x_i]} of \eqn{X}. For the weight \eqn{A_i(x)}, we use the Nadaraya-Watson weights.

Noh (2014) considered the same generalization of linear smoothers \eqn{\hat \varphi(x|p)} for frontier estimation, but with a difference choice of the weight \eqn{p}. Using the same envelopment and shape constraints as Parmeter and Racine (2013), the weight vector \eqn{p} is chosen to minimize the area under the fitted curve \eqn{\hat \varphi(x|p)}, that is \eqn{A(p) = \int_a^b\hat \varphi(x|p) dx = \sum_{i=1}^n p_i y_i \left( \int_a^b A_i(x) dx \right)}, where \eqn{[a,b]} is the true support of \eqn{X}. In practice, we integrate over the observed support \eqn{[\min_i x_i,\max_i x_i]} since the theoretic one is unknown. In what concerns the kernel weights \eqn{A_i(x)}, we use the Priestley-Chao weights 
\deqn{
A_i(x) = \left\{ \begin{array}{cc} 0 &,~i=1 \\ (x_i - x_{i-1}) K_i(x) &,~i \neq 1 \end{array} \right.,
}
where it is assumed that the pairs \eqn{(x_i,y_i)} have been ordered so that \eqn{x_1 \leq \cdots \leq x_n}. The choice of such weights is motivated by their convenience for the evaluation of the integral \eqn{\int A_i(x) dx}.
}

\value{
Returns a numeric vector with the same length as \code{x}. Returns a vector of NA if no solution has been found by the solver (GLPK). 
}

\references{
Noh, H. (2014). Frontier estimation using kernel smoothing estimators with data transformation. \emph{Journal of the Korean Statistical Society}, 43, 503-512.

Parmeter, C.F. and Racine, J.S. (2013). Smooth constrained frontier analysis in \emph{Recent Advances and Future Directions in Causality, Prediction, and Specification Analysis}, Springer-Verlag, New York, 463-488.
}

\author{
Hohsuk Noh
}

\seealso{
\code{\link{kern_smooth_bw}}.
}

\examples{
\dontrun{
data("green")
x.green <- seq(min(log(green$COST)), max(log(green$COST)), 
 length.out = 101)
options(np.tree=TRUE, crs.messages=FALSE, np.messages=FALSE)
# 1. Unconstrained 
(h.bic.green.u <- kern_smooth_bw(log(green$COST), 
 log(green$OUTPUT), method = "u", technique = "noh", 
 bw_method = "bic"))
y.ks.green.u <- kern_smooth(log(green$COST), 
 log(green$OUTPUT), x.green, h = h.bic.green.u, 
 method = "u", technique = "noh")

# 2. Monotonicity constraint
(h.bic.green.m <- kern_smooth_bw(log(green$COST),
 log(green$OUTPUT), method = "m", technique = "noh",
 bw_method = "bic"))
y.ks.green.m <- kern_smooth(log(green$COST), 
 log(green$OUTPUT), x.green, h = h.bic.green.m, 
 method = "m", technique = "noh")

# 3. Monotonicity and Concavity constraints
(h.bic.green.mc<-kern_smooth_bw(log(green$COST), log(green$OUTPUT), 
 method="mc", technique="noh", bw_method="bic"))
y.ks.green.mc<-kern_smooth(log(green$COST), 
 log(green$OUTPUT), x.green, h=h.bic.green.mc, method="mc", 
 technique="noh")

# Representation 
plot(log(OUTPUT)~log(COST), data=green, xlab="log(COST)", 
 ylab="log(OUTPUT)") 
lines(x.green, y.ks.green.u, lty=1, lwd=4, col="green")
lines(x.green, y.ks.green.m, lty=2, lwd=4, col="cyan")
lines(x.green, y.ks.green.mc, lty=3, lwd=4, col="magenta")   
legend("topleft", col=c("green","cyan","magenta"), 
lty=c(1,2,3), legend=c("unconstrained", "monotone", 
 "monotone + concave"), lwd=4, cex=0.8)
}    
}


\keyword{optimize}