\name{nuclear}
\alias{nuclear}
\docType{data}
\title{Reliability programs of nuclear reactors
}
\description{
The dataset from the US Electric Power Research Institute (EPRI) consists of 254 toughness results obtained from
non-irradiated representative steels. For each steel \eqn{i}, fracture toughness \eqn{y_i}{yi} and temperature \eqn{x_i}{xi}
were measured.
}
\usage{data(nuclear)}
\format{
  A data frame with 254 observations on the following 2 variables.
  \describe{
    \item{\code{xtab}}{Temperature.}
    \item{\code{ytab}}{Fracture toughness of each material.}
  }
}

\source{
US Electric Power Research Institute (EPRI).
}
\references{
Daouia, A., Girard, S. and Guillou, A. (2014). A Gamma-moment approach to monotonic boundary estimation. \emph{Journal of Econometrics}, 78, 727-740. 
}
\examples{
data("nuclear")
}
\keyword{datasets}
