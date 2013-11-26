\name{nuclear}
\alias{nuclear}
\docType{data}
\title{Reliability programs of nuclear reactors
}
\description{
The dataset from
the US Electric Power Research Institute (EPRI) consists of 254 toughness results obtained from
non-irradiated representative steels. For each steel \eqn{i}, fracture toughness \eqn{x_i}{xi} and temperature \eqn{y_i}{yi}
were measured.
}
\usage{data(nuclear)}
\format{
  A data frame with 254 observations on the following 2 variables.
  \describe{
    \item{\code{xtab}}{Temperature}
    \item{\code{ytab}}{fracture toughness of each material}
  }
}

\source{
US Electric Power Research Institute (EPRI)
}
\references{
Daouia, A., Girard, S. and Guillou, A. (2013). A Gamma-moment approach to monotonic boundary estimation. \emph{Journal of Econometrics}. 
http://dx.doi.org/10.1016/j.jeconom.2013.10.013 
}
\examples{
data("nuclear")
}
\keyword{datasets}
