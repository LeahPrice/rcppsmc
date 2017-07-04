#' Simple Linear Regression - Likelihood Annealing SMC
#'
#' The \code{LinReg} function provides a simple example for \pkg{RcppSMC} based on
#' estimating the parameters of a linear regression model using likelihood annealing sequential Monte Carlo.
#' 
#' @param dat		A two-column matrix or dataframe containing x and y
#' values. The data set from Williams (1959) is used as the
#' default if no data is supplied.
#' 
#' @param temperatures	Here the targets are defined by \eqn{P(y|\theta)^{\gamma_t}P(\theta)}
#' where \eqn{0=\gamma_0\le \ldots \le \gamma_T = 1} can be referred to
#' as the temperatures, \eqn{P(y|\theta)} is the likelihood and \eqn{P(\theta)}
#' is the prior. An example temperature schedule is used if one is not specified.
#' 
#' @param particles	An integer specifying the number of particles.
#' 
#' @return			A \code{list} containing the means and variances of 
#' the estimated \eqn{alpha}{alpha}, \eqn{beta}{beta} and \eqn{phi}{phi} at each iteration. The
#' effective sample size at each of the iterations and the logarithm of the estimated model
#' evidence (normalising constant of the posterior) are also returned based on several different
#' estimators.
#'
#' @details The \code{LinRegLA} function implements a simple example based on
#' estimating the parameters of a linear regression model using likelihood
#' annealing SMC. The details relevant to the default choices can be found
#' on page 49 of \url{https://eprints.qut.edu.au/101729/7/101729.pdf}.
#' Details of the original experiment can be found in Williams (1959).
#' 
#' The log evidence is -309.56 rounded to two decimal places.
#' 
#' @seealso \code{\link{LinReg}} for an implementation using data annealing SMC.
#' 
#' @references	Williams, E. (1959). Regression analysis. Wiley.
#' 
#' @examples	res <- LinRegLA(particles=1000)
#' 
#' @author		Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
#' @keywords	programming
#' @concept		static Bayesian annealing path
#' @export
LinRegLA<- function(dat, temperatures, particles=1000) {

    # if no data supplied, use default
    if (missing(dat)){
		data(radiata)
		dat <- radiata
	}
    if (missing(temperatures)) temperatures <- seq(0,1,0.05)^5

    # more eloquent tests can be added
    stopifnot(nrow(dat) > 0,
              ncol(dat) == 2,
              colnames(dat) == c("x", "y"),
			  nrow(temperatures) > 0)

    res <- LinRegLA_cpp(as.matrix(dat), as.matrix(temperatures), particles)

    invisible(res)
}