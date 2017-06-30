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
