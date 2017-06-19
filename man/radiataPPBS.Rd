\name{radiataPPBS}
\alias{radiataPPBS}
\title{Simple Linear Regression - Likelihood Annealing SMC}
\description{
  The \code{radiataPPBS} function implements a simple example based on
  estimating the parameters of a linear regression model using likelihood
  annealing SMC. The details relevant to this implementation
  can be found on page 49 of \url{https://eprints.qut.edu.au/101729/7/101729.pdf}.
  Details of the original experiment can be found in Williams (1959).
  Another implementation using data annealing SMC is available
  under \code{radiataBS}.
  
  The log evidence is -309.56 rounded to two decimal places. 
}
\usage{
  radiataPPBS(data, temperatures, particles=1000)
}
\arguments{
  \item{data}{A two-column matrix or dataframe containing x and y
    values. The data set from Williams (1959) is used as the
    default if no data is supplied.}
  \item{temperatures}{Here the targets are defined by \eqn{P(y|\theta)^{\gamma_t}P(\theta)}
  where \eqn{0=\gamma_0\le \ldots \le \gamma_T = 1} can be referred to
  as the temperatures, \eqn{P(y|\theta)} is the likelihood and \eqn{P(\theta)}
  is the prior. An example temperature schedule is used if one is not specified.}
  \item{particles}{An integer specifying the number of particles.}
}
\value{
  The \code{radiataPPBS} function returns a \code{list} containing the means and variances of 
  the estimated \eqn{alpha}{alpha}, \eqn{beta}{beta} and \eqn{phi}{phi} at each iteration. The
  effective sample size at each of the iterations and the logarithm of the estimated model
  evidence (normalising constant of the posterior) are also returned.
}
\details{
  The \code{radiataPPBS} function provides a simple example for
  \pkg{RcppSMC}.
}
\references{
  Williams, E. (1959). Regression analysis. Wiley.
}
\examples{
  res <- radiataPPBS(particles=1000)
}
\author{Adam M. Johansen, Dirk Eddelbuettel and Leah South}
\keyword{programming}
