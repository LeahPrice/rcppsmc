\name{PMMH}
\alias{PMMH}
\title{Particle marginal Metropolis-Hastings}
\description{
  The \code{PMMH} function implements particle marginal Metropolis-Hastings
  for the non-linear state space model described in Section 3.1 of
  Andrieu et al. (2010).
  
  The same prior and random walk tuning parameters are used as in Andrieu et al. (2010).
  The default data was also simulated based on the specifications in this paper,
  and the \code{simPMMH.R} script was used to simulate this data.  
  
  This toy non-linear state space model has also been described in Gordon et al. (1993)
  and Kitagawa (1996).
}
\usage{
  PMMH(data, particles=1000, iterations=10000)
}
\arguments{
  \item{data}{A vector of the observed data - the default }
  \item{particles}{An integer specifying the number of particles in the particle
  filtering estimates of the likelihood}
  \item{iterations}{An integer specifying the number of MCMC iterations}
}
\value{
  The \code{PMMH} function returns a \code{data.frame} containing the chain
  of simulated \eqn{\sigma_v} and \eqn{\sigma_w} values, as well as the
  corresponding log likelihood estimates and log prior values.
}
\details{
  The \code{PMMH} function provides an example implementation of PMMH for
  \pkg{RcppSMC}.
}
\references{
  C. Andrieu, A. Doucet, and R. Holenstein. Particle Markov chain Monte Carlo methods.
  Journal of the Royal Statistical Society: Series B (Statistical Methodology),
  72(3):269-342, 2010.
  
  N. J. Gordon, S. J. Salmond, and A. F. M. Smith. Novel approach to
  nonlinear/non-Gaussian Bayesian state estimation. IEE Proceedings-F, 
  140(2):107-113, April 1993.

  G. Kitagawa. Monte Carlo filter and smoother for non-Gaussian nonlinear
  state space models. Journal of Computational and Graphical Statistics,
  5(1):1-25, 1996.
}
\examples{
  res <- PMMH(particles=1000,iterations=10000)
}
\author{Adam M. Johansen, Dirk Eddelbuettel and Leah South}
\keyword{programming}
