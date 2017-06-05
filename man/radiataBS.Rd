\name{radiataBS}
\alias{radiataBS}
\title{Particle Filter Example}
\description{
  The \code{radiataBS} function implements a simple example based on
  estimating the parameters of a linear regression model using data
  annealing SMC. The details relevant to this implementation
  can be found on page 49 of \url{https://eprints.qut.edu.au/101729/7/101729.pdf}.
  Details of the original experiment can be found in Williams (1959).
  It would be simple to change the code to do likelihood annealing SMC.
  
  Note that there is no move step here, only reweighting and resampling.
  This example was implemented as a quick check for evidence
  estimation (the log evidence is -309.56 rounded to two decimal places).
  A MH-MCMC step will be added later to diversify the
  particles (do not expect good posterior estimates from this version). 
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
  The \code{radiataBS} function returns a \code{data.frame} containing as many rows as in
  the input data, and columns corresponding to the means and variances of 
  the estimated \eqn{alpha}{alpha}, \eqn{beta}{beta} and \eqn{phi}{phi}. The intermediate
  normalising constants, cumulative normalising constants and the ESS's are also given.
  The final cumulative normalising constant is the logarithm of the estimated evidence,
  which is useful for Bayesian model choice.
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
\author{Adam M. Johansen and Dirk Eddelbuettel}
\keyword{programming}
