\name{npbr-package}
\alias{npbr-package}
\alias{npbr}
\docType{package}
\title{Nonparametric boundary regression}
\description{
The package npbr (Daouia et al., 2017) is the first free specialized software for data edge and frontier analysis in the statistical literature.
It provides a variety of  functions for the best known and most innovative approaches to nonparametric boundary estimation. 
The selected methods are concerned with empirical, smoothed, unrestricted  as well  as constrained fits under both separate and multiple shape constraints.  
They also cover data envelopment techniques as well as robust approaches to outliers. 
The routines included in \pkg{npbr}  are  user  friendly and afford a large degree of flexibility in the estimation specifications.  
They provide smoothing parameter selection for the modern local linear and polynomial spline methods and for some promising extreme value techniques. 
Also, they seamlessly allow for Monte Carlo comparisons among the implemented estimation procedures.  This package will be very useful for statisticians and applied researchers interested 
in employing nonparametric boundary regression models. Its use is illustrated with a number of empirical applications and simulated examples.
}
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
  
where the graph of \eqn{\varphi} corresponds to the locus of the curve above which the density \eqn{f} is zero. We consider the estimation of the frontier function 
\eqn{\varphi} based on the sample \eqn{\{ (x_i,y_i),~i=1,\ldots,n\}}{{(x[i],y[i]), i=1,...,n}} in the general setting where the density \eqn{f} may have sudden jumps at the frontier, 
decay to zero or rise up to infinity as it approaches its support boundary.  

The overall objective of the present package is to provide a large variety of functions for the best known approaches to nonparametric boundary regression, 
including the vast class of methods employed in both Monte Carlo 
comparisons of  Daouia et al. (2016) and Noh (2014) as well as other promising nonparametric devices, namely the extreme-value techniques 
of Gijbels and Peng (2000), Daouia et al. (2010) and Daouia et al. (2012).   
The various functions in the \pkg{npbr} package are summarized in the table below.
We are not aware of any other existing set of statistical routines more adapted to data envelope fitting and robust frontier estimation.  
Only the classical nonsmooth FDH and DEA methods can be found in some available packages dedicated to the economic literature on measurements of  the production performance of enterprises, such as the programs 
\pkg{Benchmarking} by  Bogetoft and Otto (2011) and  \pkg{FEAR} by  Wilson (2008). Other contributions to the econometric literature on frontier analysis by Parmeter and Racine (2013) can be found at \url{http://socserv.mcmaster.ca/racinej/Gallery/Home.html}.
The package \pkg{npbr} is actually the first free specialized software for the statistical literature on nonparametric frontier analysis.  The routines included in \pkg{npbr} are  user 
friendly and highly flexible  in terms of estimation specifications.  They allow the user to locate the boundary from data 
by making use of both empirical and smooth fits as well as (un)constrained  estimates under single and multiple shape constraints.    
They also integrate smoothing parameter selection for the innovative methods based on local linear techniques, polynomial splines, extreme values and kernel smoothing,  
though the proposed selection procedures can be computationally demanding.
  
In addition, the package will be extremely useful for researchers and practitioners interested in employing nonparametric boundary regression methods. 
On one hand, such methods are very appealing because they rely on very few assumptions and benefit from their modeling flexibility, function approximation 
power and ability to detect the boundary structure of data without recourse to any \emph{a priori} parametric restrictions on the shape of the frontier and/or the distribution of noise. 
On the other hand, the package offers \R users and statisticians in this active area of research simple functions to compute the empirical mean integrated squared error, 
the empirical integrated squared bias and the empirical integrated variance of the implemented frontier estimators. 
This seamlessly allows  the interested researcher  to reproduce the Monte Carlo estimates obtained in the original articles and, perhaps most importantly,  
to easily compare the quality of any new proposal with the competitive existing methods.

 \tabular{lll}{
Function                              \tab  Description                           \tab Reference \cr
\code{\link{dea_est}}                 \tab  DEA, FDH                              \tab  Farrell (1957) \cr
                                      \tab                                        \tab  Deprins et al. (1984), \cr
                                      \tab  and linearized FDH                    \tab Hall and Park (2002) \cr
                                      \tab                                        \tab Jeong and Simar (2006)\cr            
\code{\link{loc_est}}                 \tab  Local linear fitting                  \tab  Hall et al. (1998), \cr
                                      \tab                                        \tab  Hall and Park (2004)  \cr   
\code{\link{loc_est_bw}}              \tab  Bandwidth choice                      \tab  Hall and Park (2004) \cr 
                                      \tab  for local linear fitting              \tab    \cr 
\code{\link{poly_est}}                \tab  Polynomial estimation                 \tab  Hall et al. (1998)  \cr
\code{\link{poly_degree}}             \tab Optimal polynomial                     \tab  Daouia et al. (2015)  \cr
                                      \tab  degree selection                      \tab   \cr
\code{\link{dfs_momt}}                \tab  Moment type estimation                \tab  Daouia et al. (2010),  \cr
                                      \tab                                        \tab  Dekkers et al. (1989)  \cr
\code{\link{dfs_pick}}                \tab  Pickands type estimation              \tab Daouia et al. (2010),  \cr
                                      \tab                                        \tab  Dekkers and de Haan (1989) \cr
\code{\link{rho_momt_pick}}           \tab  Conditional tail                      \tab Daouia et al. (2010),   \cr
                                      \tab index estimation                       \tab Dekkers et al. (1989),  \cr
                                      \tab                                        \tab  Dekkers and de Haan (1989) \cr
\code{\link{kopt_momt_pick}}          \tab Threshold selection for                \tab Daouia et al. (2010)  \cr
                                      \tab moment/Pickands frontiers              \tab  \cr
\code{\link{dfs_pwm}}                \tab Nonparametric frontier                 \tab Daouia et al. (2012) \cr
                                      \tab regularization                         \tab  \cr
\code{\link{loc_max}}                 \tab  Local constant estimation             \tab Gijbels and Peng (2000) \cr
\code{\link{pick_est}}                \tab  Local extreme-value estimation        \tab  Gijbels and Peng (2000) \cr
\code{\link{quad_spline_est}}         \tab  Quadratic spline fitting              \tab Daouia et al. (2015)  \cr
\code{\link{quad_spline_kn}}          \tab  Knot selection for                    \tab  Daouia et al. (2015)   \cr
                                      \tab  quadratic spline fitting              \tab    \cr
\code{\link{cub_spline_est}}          \tab  Cubic spline fitting                  \tab Daouia et al. (2015) \cr
\code{\link{cub_spline_kn}}           \tab  Knot selection for                    \tab  Daouia et al. (2015)  \cr
                                      \tab  cubic spline fitting                  \tab     \cr
\code{\link{kern_smooth}}                   \tab Nonparametric kernel                   \tab  Parmeter and Racine (2013), \cr
                                      \tab boundary regression                    \tab  Noh (2014) \cr
\code{\link{kern_smooth_bw}}               \tab Bandwidth choice for                   \tab  Parmeter and Racine (2013),    \cr
                                      \tab kernel boundary regression             \tab Noh (2014)   \cr
}

}
\author{
Abdelaati Daouia <Abdelaati.Daouia@tse-fr.eu>, Thibault Laurent <thibault.laurent@univ-tlse1.fr>, Hohsuk Noh <word5810@gmail.com>

Maintainer: Thibault Laurent <thibault.laurent@univ-tlse1.fr>
}
\references{
Daouia, A., Florens, J.-P. and Simar, L. (2010). Frontier estimation and extreme value theory. \emph{Bernoulli}, \bold{16}, 1039-1063.

Daouia, A., Florens, J.-P. and Simar, L. (2012). Regularization of Nonparametric Frontier Estimators. \emph{Journal of Econometrics}, \bold{168}, 285-299.

Daouia, A., Laurent, T. and Noh, H. (2017). npbr: A Package for Nonparametric Boundary Regression in R. \emph{Journal of Statistical Software}, \bold{79}(9), 1-43. doi:10.18637/jss.v079.i09.
	
Daouia, A., Noh, H. and Park, B.U. (2016). Data Envelope fitting with constrained polynomial splines. \emph{Journal of the Royal Statistical Society: Series B}, \bold{78}(1), 3-30. doi:10.1111/rssb.12098.

Dekkers, A.L.M. and L. de Haan (1989). On the estimation of extreme-value index and large quantiles estimation. \emph{Annals of Statistics}, \bold{17}, 1795-1832.

Dekkers, A.L.M., Einmahl, J.H.J. and L. de Haan (1989). A moment estimator for the index of an extreme-value distribution. \emph{Annals of Statistics}, \bold{17}, 1833-1855.

Deprins, D., Simar, L. and Tulkens H. (1984). Measuring labor efficiency in post offices, in: M. Marchand, P. Pestieau and H. Tulkens (Eds), The performance of Public Enterprises: Concepts and Measurements. North-Holland, Amsterdam, 243-267.

Farrell, M.J. (1957). The measurement of productive efficiency. \emph{Journal of the Royal Statistical Society, Series A}, \bold{120}, 253-281.

Gijbels, I. and Peng, L. (2000). Estimation of a support curve via order statistics. \emph{Extremes}, \bold{3}, 251-277.

Hall, P., Park, B.U. and Stern, S.E. (1998). On polynomial estimators of frontiers and boundaries. \emph{Journal of Multivariate Analysis}, \bold{66}, 71-98.

Hall, P. and Park, B.U. (2004). Bandwidth choice for local polynomial estimation of smooth boundaries. \emph{Journal of Multivariate Analysis}, \bold{91}, 240-261.

Jeong, S.-O. and Simar, L. (2006). Linearly interpolated FDH efficiency score for nonconvex frontiers. \emph{Journal of Multivariate Analysis}, \bold{97}, 2141-2161.

Noh, H. (2014). Frontier Estimation using Kernel Smoothing with Data Transformation. \emph{Journal of the Korean Statistical Society}, \bold{43}, 503-512.

Parmeter, C. and Racine, J.S. (2013). Smooth Constrained Frontier Analysis. In \emph{Recent Advances and Future Directions in Causality, Prediction, and Specification Analysis: Essays in Honor of Halbert L. White, Jr.}, Springer Verlag, (X. Chen and N.R. Swanson Eds), 463-488.
}

\keyword{nonparametric}
\keyword{optimize}


