#' Simple Linear Regression - Data Annealing SMC.
#'
#' The \code{LinReg} function provides a simple example for \pkg{RcppSMC} based on
#' estimating the parameters of a linear regression model using data annealing sequential Monte Carlo.
#' 
#' @param dat		A two-column matrix or dataframe containing x and y
#' values. The data set from Williams (1959) is used as the
#' default if no data is supplied.
#' 
#' @param particles	An integer specifying the number of particles.

#' @return			A \code{list} containing the means and variances of 
#' the estimated \eqn{alpha}{alpha}, \eqn{beta}{beta} and \eqn{phi}{phi} at each iteration. The
#' effective sample size at each of the iterations and the logarithm of the estimated model
#' evidence (normalising constant of the posterior) are also returned.
#'
#' @details The \code{LinReg} function implements a simple example based on
#' estimating the parameters of a linear regression model using data
#' annealing SMC. Details of the experiment on which the default
#' data is based can be found in Williams (1959).
#' 
#' The log evidence is -309.56 rounded to two decimal places. 
#' 
#' @seealso \code{\link{LinRegLA}} for an implementation using likelihood annealing SMC.
#' 
#' @references	Williams, E. (1959). Regression analysis. Wiley.
#' 
#' @examples	res <- LinReg(particles=1000)
#' 
#' @author		Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
#' @keywords	programming
#' @concept		static Bayesian annealing
#' @export
LinReg<- function(dat, particles=1000) {

    # if no data supplied, use default
    if (missing(dat)){
		data(radiata)
		dat <- radiata
	}

    # more eloquent tests can be added
    stopifnot(nrow(dat) > 0,
              ncol(dat) == 2,
              colnames(dat) == c("x", "y"))
	
	res <- LinReg_cpp(as.matrix(dat),particles)			 
	
    invisible(res)
}
