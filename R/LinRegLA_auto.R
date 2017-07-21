#' Simple Linear Regression - Likelihood Annealing SMC with adaptation
#'
#' The \code{LinReg_auto} function provides a simple example for \pkg{RcppSMC} based on
#' estimating the parameters of a linear regression model using likelihood annealing
#' sequential Monte Carlo. The temperature schedule and the random walk covariance
#' matrix are adapted online.
#' 
#' @param dat		A two-column matrix or dataframe containing x and y
#' values. The data set from Williams (1959) is used as the
#' default if no data is supplied.
#' 
#' @param particles	An integer specifying the number of particles.
#'
#' @param resampTol Such that resampling occurs when the ESS is less than resampTol*particles.
#' \eqn{0 \leq resampTol \leq 1}
#'
#' @param tempTol Such that the conditional ESS is maintained at tempTol*particles.
#' \eqn{0 \leq tempTol \leq 1}
#' 
#' @return			A \code{list} containing the logarithm of the estimated model
#' evidence (normalising constant of the posterior), the temperature schedule and 
#' the effective sample sizes at each temperature.
#'
#' @details The \code{LinRegLA_auto} function implements a simple example based on
#' estimating the parameters of a linear regression model using likelihood
#' annealing SMC. The details relevant to the default choices can be found
#' on page 49 of \url{https://eprints.qut.edu.au/101729/7/101729.pdf}.
#' Details of the original experiment can be found in Williams (1959).
#' 
#' The log evidence is -309.56 rounded to two decimal places.
#' 
#' @seealso \code{\link{LinRegLA}} for a non-adapted implementation.
#' 
#' @references	Williams, E. (1959). Regression analysis. Wiley.
#' 
#' @examples	res <- LinRegLA_auto(particles=1000, resampTol=0.5, tempTol = 0.9)
#' 
#' @author		Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
#' @keywords	programming
#' @concept		static Bayesian annealing path adaptive
#' @export
LinRegLA_auto<- function(dat, particles=1000, resampTol=0.5, tempTol = 0.9) {

    # if no data supplied, use default
    if (missing(dat)){
		data(radiata)
		dat <- radiata
	}
    
    # more eloquent tests can be added
    stopifnot(nrow(dat) > 0,
              ncol(dat) == 2,
              colnames(dat) == c("x", "y"))

    res <- LinRegLA_auto_cpp(as.matrix(dat), particles, resampTol, tempTol)

    invisible(res)
}