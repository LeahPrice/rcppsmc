\name{radiataBS}
\alias{radiataBS}
\title{Simple Linear Regression - Data Annealing SMC}
\description{
  The \code{radiataBS} function implements a simple example based on
  estimating the parameters of a linear regression model using data
  annealing SMC. The details relevant to this implementation
  can be found on page 49 of \url{https://eprints.qut.edu.au/101729/7/101729.pdf}.
  Details of the original experiment can be found in Williams (1959).
  Another implementation using likelihood annealing SMC is available
  under \code{radiataPPBS}.
  
  The log evidence is -309.56 rounded to two decimal places. 
}
\usage{
  radiataBS(data, particles=1000)
}
\arguments{
  \item{data}{A two-column matrix or dataframe containing x and y
    values. The data set from Williams (1959) is used as the
    default if no data is supplied.}
  \item{particles}{An integer specifying the number of particles.}
}
\value{
  The \code{radiataBS} function returns a \code{list} containing the means and variances of 
  the estimated \eqn{alpha}{alpha}, \eqn{beta}{beta} and \eqn{phi}{phi} at each iteration. The
  effective sample size at each of the iterations and the logarithm of the estimated model
  evidence (normalising constant of the posterior) are also returned.
}
\details{
  The \code{radiataBS} function provides a simple example for
  \pkg{RcppSMC}.
}
\references{
  Williams, E. (1959). Regression analysis. Wiley.
}
\examples{
  res <- radiataBS(particles=1000)
}
\author{Adam M. Johansen, Dirk Eddelbuettel and Leah F. South}
\keyword{programming}
