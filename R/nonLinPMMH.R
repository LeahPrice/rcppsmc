#' Particle marginal Metropolis-Hastings for a non-linear state space model.
#' 
#' The \code{nonLinPMMH} function implements particle marginal Metropolis-Hastings
#' for the non-linear state space model described in Section 3.1 of
#' Andrieu et al. (2010).
#' 
#' @param data			A vector of the observed data
#'
#' @param particles		An integer specifying the number of particles in the particle
#' filtering estimates of the likelihood
#' 
#' @param iterations	An integer specifying the number of MCMC iterations
#'
#' @return				A \code{data.frame} containing the chain
#' of simulated \eqn{\sigma_v} and \eqn{\sigma_w} values, as well as the
#' corresponding log likelihood estimates and log prior values.
#' 
#' @details For PMMH, The same prior and random walk tuning parameters are used as in Andrieu et al. (2010).
#' The default data was also simulated based on the specifications in this paper,
#' and the \code{simNonLinPMMH.R} script was used to simulate this data.  
#' 
#' This toy non-linear state space model has also been described in Gordon et al. (1993)
#' and Kitagawa (1996).
#'
#' @references
#' C. Andrieu, A. Doucet, and R. Holenstein. Particle Markov chain Monte Carlo methods.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology),
#' 72(3):269-342, 2010.
#' 
#' N. J. Gordon, S. J. Salmond, and A. F. M. Smith. Novel approach to
#' nonlinear/non-Gaussian Bayesian state estimation. IEE Proceedings-F, 
#' 140(2):107-113, April 1993.
#' 
#' G. Kitagawa. Monte Carlo filter and smoother for non-Gaussian nonlinear
#' state space models. Journal of Computational and Graphical Statistics,
#' 5(1):1-25, 1996.
#' 
#' @examples
#' \dontrun{
#' sim <- simNonlinPMMH(len=500)
#' res <- nonLinPMMH(sim$data,particles=1000,iterations=10000)
#' }
#' 
#' @author		Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
#' @keywords	programming
#' @concept		Bayesian PMMH PF
#' @export
nonLinPMMH<- function(data, particles=5000, iterations=10000) {

    if (missing(data)) {
         warning("data argument contained no data, using data simulated from the model.")
         data <- simNonlinPMMH(len=500)$data
    }

    res <- nonLinPMMH_cpp(as.matrix(data), particles, iterations)
	
	# Need to add some plots here - probably bivariate density plots, histogram and trace plots
	# like in Figure 4(c) of Andrieu et al. (2010). Also probably some autocorrelation plots
	# of the parameters, like in Figures 5(b) and 5(d).

    invisible(res)
}
