\name{LinRegBS}
\alias{LinRegBS}
\title{Simple Linear Regression - Data Annealing SMC}
\description{
  The \code{LinRegBS} function implements a simple example based on
  estimating the parameters of a linear regression model using data
  annealing SMC. Details of the experiment on which the default
  data is based can be found in Williams (1959).
  Another implementation using likelihood annealing SMC is available
  under \code{LinRegLABS}.
  
  The log evidence is -309.56 rounded to two decimal places. 
}
\usage{
  LinRegBS(data, particles=1000)
}
\arguments{
  \item{data}{A two-column matrix or dataframe containing x and y
    values. The data set from Williams (1959) is used as the
    default if no data is supplied.}
  \item{particles}{An integer specifying the number of particles.}
}
\value{
  The \code{LinRegBS} function returns a \code{list} containing the means and variances of 
  the estimated \eqn{alpha}{alpha}, \eqn{beta}{beta} and \eqn{phi}{phi} at each iteration. The
  effective sample size at each of the iterations and the logarithm of the estimated model
  evidence (normalising constant of the posterior) are also returned.
}
\details{
  The \code{LinRegBS} function provides a simple example for
  \pkg{RcppSMC}.
}
\references{
  Williams, E. (1959). Regression analysis. Wiley.
}
\examples{
  res <- LinRegBS(particles=1000)
}
\author{Adam M. Johansen, Dirk Eddelbuettel and Leah F. South}
\keyword{programming}
