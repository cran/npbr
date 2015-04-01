\name{post}
\Rdversion{1.1}
\alias{post}
\docType{data}
\title{French postal services}
\description{The dataset \code{post} about the cost of the delivery activity of the postal services in France was first analyzed by Cazals, Florens and Simar (2002). 
There are 4,000 post offices observed in 1994. For each post office \eqn{i}, the input \eqn{x_i}{x[i]} is the labor cost measured by the quantity of labor, 
which represents more than \eqn{80\%} of the total cost of the delivery activity. The output \eqn{y_i}{y[i]} is defined as the volume of delivered mail (in number of objects).}

\usage{data(post)}
\format{
  A data frame with 4000 observations on the following 3 variables.
  \describe{
    \item{\code{ident}}{a numeric vector.}
    \item{\code{xinput}}{a numeric vector.}
    \item{\code{yprod}}{a numeric vector.}
  }
}
\references{
Cazals, C., Florens, J.-P., Simar, L. (2002), Nonparametric frontier estimation: a robust approach, \emph{Journal of Econometrics}, 106, 1-25.
}
\examples{
data("post")
}
\keyword{datasets}
