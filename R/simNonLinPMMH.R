#' Simulates from the associated model.
#' 
#' The \code{simNonlinPMMH} function simulates data from the associated non-linear state space model.
#' 
#' @param len The length of the data sequence to simulate.
#' 
#' @return The \code{simNonlinPMMH} function returns a list containing the state and data sequences.
#' 
#' @details The \code{simNonLinPMMH} function simulates from the same model returning both
#' the state and observation vectors.
#' 
#' @rdname nonLinPMMH
#' @export
simNonlinPMMH <- function(len = 500)
{
   sim <- list()

   innovations <- rnorm(len) * sqrt(10)
   sim$state[1] <- rnorm(1) * sqrt(5)
   for (i in 2:len) {
       sim$state[i] <- 0.5 * sim$state[i-1] + 25 * sim$state[i-1] /
           (1 + sim$state[i-1]^2) + 8 * cos(1.2*i) + innovations[i]
   }
   sim$data <- sim$state^2 / 20 + rnorm(len)

   invisible(sim)
}
