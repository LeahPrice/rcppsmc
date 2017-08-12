simLinGauss <- function(len = 100, A)
{
    sim <- list()

    d <- nrow(A)

    sim$state <- matrix(rep(0,len*d),nrow=len)
    innovations <- matrix(rnorm(len*d),nrow=len)
    sim$state[1,] <- rnorm(d)
    for (i in 2:len) {
        sim$state[i,] <- A%*%sim$state[i-1,] + innovations[i,]
    }
    sim$data <- sim$state + sqrt(0.25)*matrix(rnorm(len*d),nrow=len)

    invisible(sim)
}
