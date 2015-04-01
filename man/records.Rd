\name{records}
\alias{records}
\docType{data}
\title{Annual sport records 
}
\description{
The dataset \code{records} is concerned with the yearly best men's outdoor 1500m times starting from 1966. Following Jirak, Meister and Reiss (2014), the lower boundary can be 
interpreted as the best possible time for a given year. This boundary is not believed to be shape constrained and can be estimated by any unconstrained shape nonparametric method.
}
\usage{data(records)}
\format{
  A data frame with 46 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{year.}
    \item{\code{result}}{1500m record in seconds.}
  }
}


\references{
Jirak, M.,  Meister, A.  and M. Reiss (2014), Optimal adaptive estimation in nonparametric regression with one-sided errors. \emph{Annals of Statistics}, 42, 1970--2002.
}
\examples{
data("records")
}
\keyword{datasets}
