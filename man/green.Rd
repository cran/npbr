\name{green}
\alias{green}
\docType{data}
\title{
American electric utility companies
}
\description{
The dataset consists of 123 American electric utility companies. As in the set-up of Gijbels \emph{et al.} (1999), 
we used the measurements of the variables \eqn{y_i = \log(q_i)}{yi=log(qi)} and 
\eqn{x_i = \log(c_i)}{xi=log(ci)}, where \eqn{q_i}{qi} is the production output 
of the company \eqn{i} and \eqn{c_i}{ci} is the total cost involved in the production.  
For a detailed description and analysis of these data see \emph{e.g.} Christensen and Greene (1976) and Greene (1990).
}
\usage{data(green)}
\format{
  A data frame with 123 observations on the following 2 variables.
  \describe{
    \item{\code{COST}}{a numeric vector}
    \item{\code{OUTPUT}}{a numeric vector}
  }
}

\source{
Gijbels \emph{et al.} (1999).
}
\references{
Christensen, L.R. and Greene, W.H. (1976). Economies of Scale in U.S. Electric Power Generation, \emph{Journal of Political Economy}, University of Chicago Press, 84(4), 655-76.

Gijbels, I., Mammen, E., Park, B.U. and Simar, L. (1999). On estimation of monotone and
concave frontier functions. \emph{JASA}, 94, 220-228.    

Greene, W.H. (1990). A Gamma-distributed stochastic frontier model, \emph{Journal of Econometrics}, 46(1-2), 141-163.
}
\examples{
data("green")
}
\keyword{datasets}
