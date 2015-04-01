\name{air}
\alias{air}
\docType{data}
\title{European air controllers 
}
\description{
The dataset is concerned with the assessment of the efficiency of 37 European Air Controllers. The performance of each controller can be measured by its ``distance'' 
from the upper support boundary, or equivalently, the set of the most efficient controllers. This dataset is taken from Mouchart and Simar (2002). 
Here, the activity of the controllers is described by one input (an aggregate factor of different kind of labor) 
and one output (an aggregate factor of the activity produced, based on the number of  controlled air movements, the number of  controlled flight hours, etc.). 
See also Daouia, Florens and Simar (2008). 
}
\usage{data(air)}
\format{
  A data frame with 37 observations on the following 2 variables.
  \describe{
    \item{\code{xtab}}{an input.}
    \item{\code{ytab}}{an output.}
  }
}


\references{
Daouia, A., Florens, J.-P. and Simar, L. (2008). Functional Convergence of Quantile-type Frontiers with Application to Parametric Approximations. \emph{Journal of Statistical Planning and Inference}, 138, 708-725.

Mouchart, M. and L. Simar (2002). Efficiency analysis of Air Controllers:  first insights, Consulting report 0202, Institut de Statistique, UCL, Belgium.
}
\examples{
data("air")
}
\keyword{datasets}
